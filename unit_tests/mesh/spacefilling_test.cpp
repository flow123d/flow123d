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

struct PointHilbert {
    arma::vec2 coords;
    double hilbertIndex;
    uint originalIndex;
};

class PointsManager {
private:
    std::vector<PointHilbert> points;
    std::vector<double> sizes;
    arma::vec2 shift;
    double scalar;
    static bool comparePoints(PointHilbert& first, PointHilbert& second) {
        return first.hilbertIndex < second.hilbertIndex;
    }
public:
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
    std::vector<PointHilbert>& getPoints() {
        return points;
    }
    std::vector<double>& getSizes() {
        return sizes;
    }
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
    void setPoints(const Armor::Array<double>& data) {
        points.resize(data.size());
        arma::vec3 tmpVec;
        PointHilbert p;
        for (uint i = 0; i < data.size(); ++i) {
            tmpVec = data.vec<3>(i);
            p.coords[0] = tmpVec[0];
            p.coords[1] = tmpVec[1];
            p.originalIndex = i;
            points[i] = p;
        }
    }
    void setBoundingBox(double min1, double min2, double max1, double max2) {
        shift = arma::vec2{min1, min2};
        arma::vec2 tmp{max1, max2};
        tmp -= shift;
        scalar = tmp[0] > tmp[1] ? tmp[0] :tmp[1];
    }
    void calculateHibert(double eps) {
        for (auto& point : points) {
            point.hilbertIndex = hilbertIndex(point.coords[0], point.coords[1], eps);
        }
    }
    void calculateHibert() {
        for (uint i = 0; i < points.size(); ++i) {
            points[i].hilbertIndex = hilbertIndex(points[i].coords[0], points[i].coords[1], sizes[points[i].originalIndex] * sizes[points[i].originalIndex]);
        }
    }
    void sortByHilbert() {
        std::sort(points.begin(), points.end(), comparePoints);
    }
    void normalize() {
        for (auto& p : points) {
            p.coords = (p.coords - shift) / scalar;
        }
    }
    void fillGroupIndices(std::vector<double>& indices) const {
        indices.resize(points.size());
        std::array<uint, 10> groupIndices = {3, 8, 9, 4, 2, 6, 7, 0, 1, 5};
        uint groupIndex = 0;
        for (uint i = 0; i < indices.size(); ++i) {
            if (!(i % GROUP_SIZE)) {
                ++groupIndex;
            }
            indices[points[i].originalIndex] = groupIndices[groupIndex % 10];
        }
    }
    std::vector<uint> getForwardPermutation() const {
        std::vector<uint> forwardPermutation(points.size());
        for (uint i = 0; i < points.size(); ++i) {
            forwardPermutation[points[i].originalIndex] = i;
        }
        return forwardPermutation;
    }
    std::vector<uint> getBackwardPermutation() const {
        std::vector<uint> backwardPermutation(points.size());
        for (uint i = 0; i < points.size(); ++i) {
            backwardPermutation[i] = points[i].originalIndex;
        }
        return backwardPermutation;
    }
};

TEST(Spacefilling, get_centers) {
    Profiler::initialize();

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

//     std::string mesh_in_string = "{mesh_file=\"mesh/square_uniform.msh\"}";
//     std::string mesh_in_string = "{mesh_file=\"mesh/square_uniform_2.msh\"}";
//     std::string mesh_in_string = "{mesh_file=\"mesh/square_refined.msh\"}";
//     std::string mesh_in_string = "{mesh_file=\"mesh/lshape_refined.msh\"}";
    std::string mesh_in_string = "{mesh_file=\"mesh/lshape_refined_2.msh\"}";
    
    Mesh * mesh = mesh_full_constructor(mesh_in_string);
    
    auto reader = reader_constructor(mesh_in_string);
    reader->read_physical_names(mesh);
    reader->read_raw_mesh(mesh);
    
    std::vector<arma::vec3> centers;
    std::vector<double> sizes;
    
    
//     START_TIMER("get_centers");
        for (ElementAccessor<3> elm : mesh->elements_range()) {
            centers.push_back(elm.centre());
            sizes.push_back(std::min({
                arma::norm(*elm.node(0) - *elm.node(1)),
                arma::norm(*elm.node(1) - *elm.node(2)),
                arma::norm(*elm.node(2) - *elm.node(0))
            }));
        }
//     END_TIMER("get_centers");

    double checksum1 = 0;

    START_TIMER("calculation_before_sort");
    for (ElementAccessor<3> elm : mesh->elements_range()) {
        auto n0 = elm.node(0);
        auto n1 = elm.node(1);
        auto n2 = elm.node(2);
        arma::mat::fixed<2, 2> M;
        arma::vec3 tmp;
        tmp = *n1 - *n0;
        M.col(0) = arma::vec2{tmp[0], tmp[1]};
        tmp = *n2 - *n0;
        M.col(1) = arma::vec2{tmp[0], tmp[1]};
        double detM = arma::det(M);
        double jac = std::abs(detM) / 2;
        arma::mat::fixed<3, 3> phi = {{1, -1, -1}, {0, 1, 0}, {0, 0, 1}};
        arma::mat::fixed<3, 2> grad_phi = {{-1, -1}, {1, 0}, {0, 1}};
        arma::mat::fixed<3, 3> A_local = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (uint i = 0; i < 3; ++i) {
            for (uint j = 0; j < 3; ++j) {
                A_local(i, j) += arma::dot(M * grad_phi[i], M * grad_phi[j]) * jac;
                checksum1 += A_local(i, j);
            }
        }
    }
    END_TIMER("calculation_before_sort");
    
    std::cout << "checksum 1: " << checksum1 << '\n';
    
    Armor::Array<double>& nodes = mesh->nodes_;
    std::vector<Element>& elements = mesh->element_vec_;
    PointsManager pm2;
    
//     uint j = 0;
//     for (ElementAccessor<3> elm : mesh->elements_range()) {
//         if (j > 6262 && j < 6268) {
//             std::cout << j << ": " << pm2.hilbertIndex(elm.centre()[0] / 2, elm.centre()[1] / 2, 0.0000000000001) << '\n';
//         }
//         ++j;
//     }
    
    std::cout << "n_nodes: " << mesh->n_nodes() << '\n';
    std::vector<arma::vec3> nodes_backup(mesh->n_nodes());
    for (uint i = 0; i < nodes_backup.size(); ++i) {
        nodes_backup[i] = nodes.vec<3>(i);
    }
    BoundingBox bb(nodes_backup);

    pm2.setPoints(nodes_backup);
    pm2.setBoundingBox(bb.min()[0], bb.min()[1], bb.max()[0], bb.max()[1]);
    std::cout << pm2.getPoints()[50].coords[0] << ' ' << pm2.getPoints()[50].coords[1] << '\n';
    pm2.normalize();
    std::cout << pm2.getPoints()[50].coords[0] << ' ' << pm2.getPoints()[50].coords[1] << '\n';
    pm2.calculateHibert(0.0000000000001);
    pm2.sortByHilbert();
    std::vector<uint> newNodeIndices = pm2.getForwardPermutation();
//     for (uint i = 0; i < elements.size(); ++i) {
//         elements[newNodeIndices[i]] = elements_backup[i];
//     }
    for (uint i = 0; i < nodes_backup.size(); ++i) {
        nodes.set(newNodeIndices[i]) = nodes_backup[i];
    }
    
//     for (uint i = 0; i < mesh->n_nodes(); ++i) {
//         if (!(i % 10)) {
//             arma::vec3 tmp = nodes.vec<3>(i);
//             std::cout << i << ": " << pm2.hilbertIndex(tmp[0] / 2, tmp[1] / 2, 0.0000000000001) << '\n';
//         }
//     }
    
    vector<Element> elements_backup = elements;
    PointsManager pm;
    pm.setPoints(centers);
    pm.getSizes() = sizes;
    pm.setBoundingBox(bb.min()[0], bb.min()[1], bb.max()[0], bb.max()[1]);
    pm.normalize();
    pm.calculateHibert();
//     std::cout << "before sort: " << pm.hilbertIndex(pm.getPoints()[6264].coords[0], pm.getPoints()[6264].coords[1], 0.0000000000001) << '\n';
    pm.sortByHilbert();
//     std::cout << "after sort: " << pm.hilbertIndex(pm.getPoints()[8855].coords[0], pm.getPoints()[8855].coords[1], 0.0000000000001) << '\n';
    std::vector<uint> newElementIndices = pm.getForwardPermutation();
//     std::cout << elements.size() << '\n';
    for (uint i = 0; i < centers.size(); ++i) {
        elements[newElementIndices[i]] = elements_backup[i];
        elements[i].nodes_[0] = newNodeIndices[elements[i].nodes_[0]];
        elements[i].nodes_[1] = newNodeIndices[elements[i].nodes_[1]];
        elements[i].nodes_[2] = newNodeIndices[elements[i].nodes_[2]];
//         if (newElementIndices[i] == 8855) {
//             std::cout << "HERE: " << i << '\n';
//         }
    }
    
//     uint i = 0;
//     for (ElementAccessor<3> elm : mesh->elements_range()) {
//         if (i > 8852 && i < 8858) {
//             std::cout << i << ": " << pm2.hilbertIndex(elm.centre()[0] / 2, elm.centre()[1] / 2, 0.0000000000001) << '\n';
//         }
//         if (!(i % 100)) {
//             std::cout << i << ": " << pm2.hilbertIndex(elm.centre()[0] / 2, elm.centre()[1] / 2, 0.0000000000001) << '\n';
//         }
//         ++i;
//     }
    
    vector<Element> elements_backup_2 = elements;
    PointsManager pm3;
    pm3.setPoints(centers);
    pm3.getSizes() = sizes;
    pm3.setBoundingBox(bb.min()[0], bb.min()[1], bb.max()[0], bb.max()[1]);
    pm3.normalize();
    pm3.calculateHibert();
    pm3.sortByHilbert();
    
//     const auto pb = pm3.getPoints();
//     for (uint i = 0; i < pb.size(); ++i) {
//         if (!(i % 100)) {
//             std::cout << i << ": " << pb[i].hilbertIndex << '\n';
//         }
//     }
    
    double checksum2 = 0;
    
    START_TIMER("calculation_after_sort");
    for (ElementAccessor<3> elm : mesh->elements_range()) {
        auto n0 = elm.node(0);
        auto n1 = elm.node(1);
        auto n2 = elm.node(2);
        arma::mat::fixed<2, 2> M;
        arma::vec3 tmp;
        tmp = *n1 - *n0;
        M.col(0) = arma::vec2{tmp[0], tmp[1]};
        tmp = *n2 - *n0;
        M.col(1) = arma::vec2{tmp[0], tmp[1]};
        double detM = arma::det(M);
        double jac = std::abs(detM) / 2;
        arma::mat::fixed<3, 3> phi = {{1, -1, -1}, {0, 1, 0}, {0, 0, 1}};
        arma::mat::fixed<3, 2> grad_phi = {{-1, -1}, {1, 0}, {0, 1}};
        arma::mat::fixed<3, 3> A_local = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (uint i = 0; i < 3; ++i) {
            for (uint j = 0; j < 3; ++j) {
                A_local(i, j) += arma::dot(M * grad_phi[i], M * grad_phi[j]) * jac;
                checksum2 += A_local(i, j);
            }
        }
    }
    END_TIMER("calculation_after_sort");
    
    std::cout << "checksum 2: " << checksum2 << '\n';
    
    std::vector<double> groupIndices;
    
    pm.fillGroupIndices(groupIndices);
    
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
        data_cache.store_value(el.idx(), &groupIndices[el.idx()]);
    }
    output->write_time_frame();

    delete mesh;
    
    Profiler::instance()->output(cout);
    Profiler::uninitialize();
}

