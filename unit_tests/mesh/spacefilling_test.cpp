/*
 * field_speed_test.cpp
 *
 *  Created on: Feb 3, 2013
 *      Author: jb
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>

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

#define GROUP_SIZE 64
#define EPS 0.0001

struct PointHilbert {
    arma::vec2 coords;
    double hilbertIndex;
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
};

static const int loop_call_count = 100000;

TEST(Spacefilling, get_centers) {
    Profiler::initialize();
    
    std::vector<arma::vec3> points;

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    std::string mesh_in_string = "{mesh_file=\"mesh/square_uniform.msh\"}";
//     std::string mesh_in_string = "{mesh_file=\"mesh/square_refined.msh\"}";
//     std::string mesh_in_string = "{mesh_file=\"mesh/lshape_refined.msh\"}";
    
    Mesh * mesh = mesh_constructor(mesh_in_string);
    
    auto reader = reader_constructor(mesh_in_string);
    reader->read_physical_names(mesh);
    reader->read_raw_mesh(mesh);
    
    START_TIMER("get_centers");
//     for (int i = 0; i < loop_call_count; ++i) {
        for (ElementAccessor<3> elm : mesh->elements_range()) {
            points.push_back(elm.centre());
        }
//     }
    END_TIMER("get_centers");
    
    points.front().print();
    
    PointsManager pm;
    pm.setPoints(points);
    
    START_TIMER("calculate_hibert");
    pm.calculateHibert(EPS);
    END_TIMER("calculate_hibert");
    
    START_TIMER("sort_by_hibert");
    pm.sortByHilbert();
    END_TIMER("sort_by_hibert");
    
    pm.getPoints().front().coords.print();

    delete mesh;
    
    Profiler::instance()->output(cout);
    Profiler::uninitialize();
}

