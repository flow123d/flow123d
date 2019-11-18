/*
 * field_speed_test.cpp
 *
 *  Created on: Feb 3, 2013
 *      Author: jb
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest_mpi.hh>

#include "system/global_defs.h"

#include <mesh_constructor.hh>
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"
#include "io/msh_gmshreader.h"

#include <iostream>

#include "fields/generic_field.hh"
#include "io/output_mesh.hh"
#include "io/output_time.hh"


#define GROUP_SIZE 64
#define EPS 0.0000000001

struct PointHilbert {
    arma::vec2 coords;
    double hilbertIndex;
    uint originalIndex;
};

class PointsManager {
private:
    std::vector<PointHilbert> points;
//     std::array<arma::vec3, 6> colors;
    double hilbertIndex(double x, double y, double eps) {
        if (eps > 1) {
            return 0;
        } else {
            if (x < 0.5) {
                if (y < 0.5) {
                    return hilbertIndex(2 * y, 2 * x, 4 * eps) / 4;
                } else {
                    return (1 + hilbertIndex(2 * x, 2 * y - 1, 4 * eps)) / 4;
                }
            } else {
                if (y >= 0.5) {
                    return (2 + hilbertIndex(2 * x - 1, 2 * y - 1, 4 * eps)) / 4;
                } else {
                    return (3 + hilbertIndex(1 - 2 * y, 2 - 2 * x, 4 * eps)) / 4;
                }
            }
        }
    }
    static bool comparePoints(PointHilbert& first, PointHilbert& second) {
        return first.hilbertIndex < second.hilbertIndex;
    }
public:
//     PointsManager() {
//         colors[0] = {0, 0, 1};
//         colors[1] = {0, 1, 0};
//         colors[2] = {0, 1, 1};
//         colors[3] = {1, 0, 0};
//         colors[4] = {1, 0, 1};
//         colors[5] = {1, 1, 0};
//     }
    std::vector<PointHilbert>& getPoints() {
        return points;
    }
//     std::array<arma::vec3, 6>& getColors() {
//         return colors;
//     }
//     void generate(uint size) {
//         points.resize(size);
//         arma::vec2 vec;
//         for (auto& point : points) {
//             vec.randu();
//             point.coords = vec;
//         }
//     }
    void setPoints(const std::vector<arma::vec3>& data) {
        points.resize(data.size());
        PointHilbert p;
        for (uint i = 0; i < data.size(); ++i) {
            p.coords[0] = data[i][0];
            p.coords[1] = data[i][1];
            p.originalIndex = i;
            points[i] = p;
        }
    }
    void calculateHibert(double eps) {
        for (auto& point : points) {
            point.hilbertIndex = hilbertIndex(point.coords[0], point.coords[1], eps);
        }
    }
    void sortByHilbert() {
        std::sort(points.begin(), points.end(), comparePoints);
    }
    void fillIndices(std::vector<double>& indices) const {
        indices.resize(points.size());
        std::array<uint, 10> groupIndices = {3, 8, 9, 4, 2, 6, 7, 0, 1, 5};
        uint groupIndex = 0;
        for (uint i = 0; i < indices.size(); ++i) {
            if (!(i % 64)) {
                ++groupIndex;
            }
            indices[points[i].originalIndex] = groupIndices[groupIndex % 10];
//             indices[points[i].originalIndex] = (groupIndex % 10);
        }
    }
};

TEST(Spacefilling, get_centers) {
    Profiler::initialize();
    
    std::vector<arma::vec3> points;

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

//     std::string mesh_in_string = "{mesh_file=\"mesh/square_uniform.msh\"}";
//     std::string mesh_in_string = "{mesh_file=\"mesh/square_refined.msh\"}";
    std::string mesh_in_string = "{mesh_file=\"mesh/lshape_refined.msh\"}";
    
    Mesh * mesh = mesh_full_constructor(mesh_in_string);
    
    auto reader = reader_constructor(mesh_in_string);
    reader->read_physical_names(mesh);
    reader->read_raw_mesh(mesh);
    
    START_TIMER("get_centers");
        for (ElementAccessor<3> elm : mesh->elements_range()) {
            points.push_back(elm.centre());
        }
    END_TIMER("get_centers");
    
//     points.front().print();
    
    PointsManager pm;
    pm.setPoints(points);
    
    START_TIMER("calculate_hibert");
    pm.calculateHibert(EPS);
    END_TIMER("calculate_hibert");
    
    START_TIMER("sort_by_hibert");
    pm.sortByHilbert();
    END_TIMER("sort_by_hibert");
    
//     pm.getPoints().front().coords.print();
    
    std::vector<double> hilbert_id;
    
    pm.fillIndices(hilbert_id);
    
//     for (uint i = 0; i < 20; ++i) {
//         std::cout << hilbert_id[i] << '\n';
//     }
    
    const string test_output_time_input = R"JSON(
    {
        format = {
            TYPE = "gmsh",
            variant = "ascii"
        }
    }
    )JSON";

    auto in_rec = Input::ReaderToStorage(test_output_time_input,
    const_cast<Input::Type::Record &>(OutputTime::get_input_type()),
    Input::FileFormat::format_JSON).get_root_interface<Input::Record>();
    auto output = OutputTime::create_output_stream("dummy_equation", in_rec, "s");
    std::shared_ptr<OutputMeshBase> output_mesh = std::make_shared<OutputMeshDiscontinuous>(*mesh);
    output_mesh->create_sub_mesh();
    output->set_output_data_caches(output_mesh);
    output->update_time(0.0);
    auto data_cache = output->prepare_compute_data<double>("patch_id", OutputTime::ELEM_DATA, 1, 1);

    for(auto el : mesh->elements_range()) {
        data_cache.store_value(el.idx(), &hilbert_id[el.idx()]);
    }
    output->write_time_frame();

    delete mesh;
    
    Profiler::instance()->output(cout);
    Profiler::uninitialize();
}

