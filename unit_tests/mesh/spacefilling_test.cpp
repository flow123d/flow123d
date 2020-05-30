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
#include <algorithm>

#include "fields/generic_field.hh"
#include "io/output_mesh.hh"
#include "io/output_time.hh"


// #define GROUP_SIZE 64

double petscCalculation(Mesh * mesh) {
    PetscInt dofs[4];
    Mat global_matrix;
    int n_global_dofs = mesh->n_nodes();
    PetscInt *nnz = new PetscInt[n_global_dofs];
    // determine preallocation size
    for (int i = 0; i < n_global_dofs; ++i) {
        nnz[i] = 0;
    }
    for (ElementAccessor<3> elm : mesh->elements_range()) {
        for(uint i = 0; i <= elm.dim(); ++i) {
            uint i_dof = elm.node(i).idx();
            nnz[i_dof] += 1;
        }
    }
//     for (int i = 0; i < n_global_dofs; ++i) {
//         nnz[i] += 1;
//     }
//     MatCreateSeqAIJ(PETSC_COMM_SELF, n_global_dofs, n_global_dofs, 0, nnz, &global_matrix);
    MatCreateSeqAIJ(PETSC_COMM_SELF, n_global_dofs, n_global_dofs, 30, NULL, &global_matrix);
//     MatSetOption(global_matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatAssemblyBegin(global_matrix, MAT_FINAL_ASSEMBLY);
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
        arma::mat::fixed<4, 4> local_matrix = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
        for (uint i = 0; i < 4; ++i) {
            for (uint j = 0; j < 4; ++j) {
                local_matrix(i, j) += arma::dot(M * grad_phi[i], M * grad_phi[j]) * jac;
            }
        }
        for(uint i = 0; i <= elm.dim(); ++i) {
            dofs[i] = elm.node(i).idx();
        }
        if (MatSetValues(global_matrix, 4, dofs, 4, dofs, &local_matrix[0], ADD_VALUES) != 0) {
            break;
        }
    }
    MatAssemblyEnd(global_matrix, MAT_FINAL_ASSEMBLY);
    PetscScalar result;
    int idx = n_global_dofs - 1;
    MatGetValues(global_matrix, 1, &idx, 1, &idx, &result);
    return result;
}

Mat getGlobalMatrix(Mesh * mesh) {
    PetscInt dofs[4];
    Mat global_matrix;
    int n_global_dofs = mesh->n_nodes();
    PetscInt *nnz = new PetscInt[n_global_dofs];
    // determine preallocation size
    for (int i = 0; i < n_global_dofs; ++i) {
        nnz[i] = 0;
    }
    for (ElementAccessor<3> elm : mesh->elements_range()) {
        for(uint i = 0; i <= elm.dim(); ++i) {
            uint i_dof = elm.node(i).idx();
            nnz[i_dof] += 1;
        }
    }
//     MatCreateSeqAIJ(PETSC_COMM_SELF, n_global_dofs, n_global_dofs, 0, nnz, &global_matrix);
    MatCreateSeqAIJ(PETSC_COMM_SELF, n_global_dofs, n_global_dofs, 30, NULL, &global_matrix);
//     MatSetOption(global_matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    MatAssemblyBegin(global_matrix, MAT_FINAL_ASSEMBLY);
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
        arma::mat::fixed<4, 4> local_matrix = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
        for (uint i = 0; i < 4; ++i) {
            for (uint j = 0; j < 4; ++j) {
                local_matrix(i, j) += arma::dot(M * grad_phi[i], M * grad_phi[j]) * jac;
            }
        }
        for(uint i = 0; i <= elm.dim(); ++i) {
            dofs[i] = elm.node(i).idx();
        }
        if (MatSetValues(global_matrix, 4, dofs, 4, dofs, &local_matrix[0], ADD_VALUES) != 0) {
            break;
        }
    }
    MatAssemblyEnd(global_matrix, MAT_FINAL_ASSEMBLY);
    return global_matrix;
}

double performMultiplication(Mat& global_matrix, PetscInt n_global_dofs) {
    Vec x, y;
    VecCreateSeq(PETSC_COMM_SELF, n_global_dofs, &x);
    VecDuplicate(x, &y);
    for (uint i = 0; i < 200; ++i) {
        MatMult(global_matrix, x, y);
    }
    PetscInt ix = 0;
    PetscScalar result;
    VecGetValues(y, 1, &ix, &result);
    return result;
}

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

TEST(Spacefilling, space_filling) {
    Profiler::initialize();

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    const std::string meshName = "square_uniform_3D";
//     const std::string meshName = "square_refined_3D";
//     const std::string meshName = "lshape_uniform_3D";
//     const std::string meshName = "lshape_refined_3D";
    
//     const std::string meshName = "square_uniform_2D";
//     const std::string meshName = "square_refined_2D";
//     const std::string meshName = "lshape_uniform_2D";
//     const std::string meshName = "lshape_refined_2D";
    
//     const std::string meshName = "square_uniform_3D_small";
//     const std::string meshName = "square_refined_3D_small";
//     const std::string meshName = "lshape_uniform_3D_small";
//     const std::string meshName = "lshape_refined_3D_small";
    
//     const std::string meshName = "square_uniform_2D_small";
//     const std::string meshName = "square_refined_2D_small";
//     const std::string meshName = "lshape_uniform_2D_small";
//     const std::string meshName = "lshape_refined_2D_small";
//  --------------------------------------------------------
//     const std::string testName = "first";
//     const std::string testName = "mean";
//     const std::string testName = "zcurve";
    const std::string testName = "hilbert";
    
//     const std::string testName = "first_petsc";
//     const std::string testName = "mean_petsc";
//     const std::string testName = "zcurve_petsc";
//     const std::string testName = "hilbert_petsc";
    
//     const std::string testName = "first_petsc_matmult";
//     const std::string testName = "mean_petsc_matmult";
//     const std::string testName = "zcurve_petsc_matmult";
//     const std::string testName = "hilbert_petsc_matmult";
    
    const std::string mesh_in_string = "{mesh_file=\"mesh/" + meshName + ".msh\"}";
    
    Mesh * mesh = mesh_full_constructor(mesh_in_string);
    
    std::cout << "copying mesh" << '\n';
    Mesh * meshBackup = mesh_full_constructor(mesh_in_string);
    
//     Mat global_matrix_before_sort = getGlobalMatrix(mesh);
    
//     MeshOptimizer<2> mo(*mesh);
    MeshOptimizer<3> mo(*mesh);
    std::cout << "calculating sizes" << '\n';
    mo.calculateSizes();
    std::cout << "calculating node curve values" << '\n';
    mo.calculateNodeCurveValuesAsHilbert();
//     mo.calculateNodeCurveValuesAsZCurve();
//     mo.calculateNodeCurveValuesAsMeanOfCoords();
//     mo.calculateNodeCurveValuesAsFirstCoord();
    std::cout << "calculating element curve values" << '\n';
    mo.calculateElementCurveValuesAsHilbertOfCenters();
//     mo.calculateElementCurveValuesAsZCurveOfCenters();
//     mo.calculateElementCurveValuesAsMeanOfCoords();
//     mo.calculateElementCurveValuesAsFirstCoord();
//     mo.calculateElementCurveValuesAsMeanOfNodes();
    std::cout << "sorting nodes" << '\n';
    mo.sortNodes();
    std::cout << "sorting elements" << '\n';
    mo.sortElements();
    
    std::cout << "performing calculation" << '\n';
    START_TIMER("calculation_after_sort");
//     double result1 = calculation2DBeforeSort(mesh);
    double result1 = calculation3DBeforeSort(mesh);
//     double result1 = petscCalculation(mesh);
//     double result1 = performMultiplication(global_matrix_before_sort, mesh->n_nodes());
    END_TIMER("calculation_after_sort");
    
    std::cout << "result 1: " << result1 << '\n';
    
//     Mat global_matrix_after_sort = getGlobalMatrix(mesh);
    std::cout << "performing calculation" << '\n';
    START_TIMER("calculation_before_sort");
//     double result2 = calculation2DAfterSort(mesh);
    double result2 = calculation3DAfterSort(meshBackup);
//     double result2 = petscCalculation(mesh);
//     double result2 = performMultiplication(global_matrix_after_sort, mesh->n_nodes());
    END_TIMER("calculation_before_sort");
    
    std::cout << "result 2: " << result2 << '\n';
    
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

