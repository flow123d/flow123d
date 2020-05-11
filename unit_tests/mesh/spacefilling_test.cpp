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
#include "mesh/mesh_optimizer.hh"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"
#include "io/msh_gmshreader.h"

#include <iostream>

#include "fields/generic_field.hh"
#include "io/output_mesh.hh"
#include "io/output_time.hh"


// #define GROUP_SIZE 64

double calculation2DBeforeSort(Mesh * mesh) {
    double checksum = 0;
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
                    checksum += std::abs(A_local(i, j));
                }
            }
        }
    return checksum;  
}

double calculation2DAfterSort(Mesh * mesh) {
    double checksum = 0;
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
                    checksum += std::abs(A_local(i, j));
                }
            }
        }
    return checksum;  
}

double calculation3DBeforeSort(Mesh * mesh) {
    double checksum = 0;
    for (ElementAccessor<3> elm : mesh->elements_range()) {
        auto n0 = elm.node(0);
        auto n1 = elm.node(1);
        auto n2 = elm.node(2);
        auto n3 = elm.node(3);
        arma::mat::fixed<3, 3> M;
        M.col(0) = *n1 - *n0;
        M.col(1) = *n2 - *n0;
        M.col(2) = *n3 - *n0;
        double detM = arma::det(M);
        double jac = std::abs(detM) / 6;
        arma::mat::fixed<4, 4> phi_coefs = {{1, -1, -1,-1}, {0, 1, 0,0}, {0, 0, 1,0},{0, 0, 0, 1}};
        arma::mat::fixed<4, 3> grad_phi = {{-1, -1, -1}, {1, 0,0}, {0, 1,0}, {0, 0, 1}};
        arma::mat::fixed<4, 4> A_local = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
        for (uint i = 0; i < 4; ++i) {
            for (uint j = 0; j < 4; ++j) {
                A_local(i, j) += arma::dot(M * grad_phi[i], M * grad_phi[j]) * jac;
                checksum += std::abs(A_local(i, j));
            }
        }
    }
    return checksum;
}

double calculation3DAfterSort(Mesh * mesh) {
    double checksum = 0;
    for (ElementAccessor<3> elm : mesh->elements_range()) {
        auto n0 = elm.node(0);
        auto n1 = elm.node(1);
        auto n2 = elm.node(2);
        auto n3 = elm.node(3);
        arma::mat::fixed<3, 3> M;
        M.col(0) = *n1 - *n0;
        M.col(1) = *n2 - *n0;
        M.col(2) = *n3 - *n0;
        double detM = arma::det(M);
        double jac = std::abs(detM) / 6;
        arma::mat::fixed<4, 4> phi_coefs = {{1, -1, -1,-1}, {0, 1, 0,0}, {0, 0, 1,0},{0, 0, 0, 1}};
        arma::mat::fixed<4, 3> grad_phi = {{-1, -1, -1}, {1, 0,0}, {0, 1,0}, {0, 0, 1}};
        arma::mat::fixed<4, 4> A_local = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
        for (uint i = 0; i < 4; ++i) {
            for (uint j = 0; j < 4; ++j) {
                A_local(i, j) += arma::dot(M * grad_phi[i], M * grad_phi[j]) * jac;
                checksum += std::abs(A_local(i, j));
            }
        }
    }
    return checksum;
}

// void printMesh(Mesh& mesh) {
//     arma::vec3 sum = {0, 0, 0};
//     uint j = 0;
//     for (ElementAccessor<3> elm : mesh.elements_range()) {
//         const Element& el = *elm.element();
//         std::cout << "Element " << j << ":\n";
//         for(uint i = 0; i < el.n_nodes(); ++i) {
//             arma::vec3 tmp = mesh.nodes_.vec<3>(el.nodes_[i]);
//             sum += tmp;
//             std::cout << "    " << tmp[0] << ' ' << tmp[1] << ' ' << tmp[2] << '\n';
//         }
//         ++j;
//     }
//     std::cout << "SUM = " << sum[0] << ' ' << sum[1] << ' ' << sum[2] << '\n';
// }

TEST(Spacefilling, space_filling) {
    Profiler::initialize();

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

//     const std::string testName = "square_uniform;
//     const std::string testName = "square_uniform_2";
//     const std::string testName = "square_refined";
//     const std::string testName = "lshape_refined";
    const std::string meshName = "lshape_refined_2_cube";
    const std::string testName = "ne_hilb";
    const std::string mesh_in_string = "{mesh_file=\"mesh/" + meshName + ".msh\"}";
    
    Mesh * mesh = mesh_full_constructor(mesh_in_string);
    
//     auto reader = reader_constructor(mesh_in_string);
//     reader->read_physical_names(mesh);
//     reader->read_raw_mesh(mesh);
    
//     printMesh(*mesh);

    START_TIMER("calculation_before_sort");
//     double checksum1 = calculation2DBeforeSort(mesh);
    double checksum1 = calculation3DBeforeSort(mesh);
    END_TIMER("calculation_before_sort");
    
    std::cout << "checksum 1: " << checksum1 << '\n';
    
    MeshOptimizer<3> mo(*mesh);
    std::cout << "reading nodes" << '\n';
    mo.readNodesFromMesh();
    std::cout << "reading elemenst" << '\n';
    mo.readElementsFromMesh();
    std::cout << "calculating sizes" << '\n';
    mo.calculateSizes();
    std::cout << "calculating node curve values" << '\n';
//     mo.calculateNodeCurveValuesAsHilbert();
    mo.calculateNodeCurveValuesAsMeanOfCoords();
//     mo.calculateNodeCurveValuesAsFirstCoord();
    std::cout << "calculating element curve values" << '\n';
    mo.calculateElementCurveValuesAsMeanOfNodes();
    std::cout << "sorting nodes" << '\n';
    mo.sortNodes();
    std::cout << "sorting elements" << '\n';
    mo.sortElements();
    std::cout << "writing nodes to mesh" << '\n';
    mo.writeNodesToMesh();
    std::cout << "writing elements to mesh" << '\n';
    mo.writeElementsToMesh();
    
    START_TIMER("calculation_after_sort");
//     double checksum2 = calculation2DAfterSort(mesh);
    double checksum2 = calculation3DAfterSort(mesh);
    END_TIMER("calculation_after_sort");
    
    std::cout << "checksum 2: " << checksum2 << '\n';
    
//     printMesh(*mesh);
    
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

//     std::vector<double> groupIndices;

//     pm.fillGroupIndices(groupIndices);

//     for(auto el : mesh->elements_range()) {
//         data_cache.store_value(el.idx(), &groupIndices[el.idx()]);
//     }
//     output->write_time_frame();

    delete mesh;
    
    std::time_t unixTime = std::time(nullptr);
    std::stringstream tmpStream;
    tmpStream << unixTime;
    std::string unixTimeString = tmpStream.str();
    std::ofstream jsonResult("../../../results/" + meshName + '_' + testName + '/' + unixTimeString + ".json");
    
    Profiler::instance()->output(jsonResult);
    Profiler::uninitialize();
}

