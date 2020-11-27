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
#include <array>
#include <sstream>

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
    MatCreateSeqAIJ(PETSC_COMM_SELF, n_global_dofs, n_global_dofs, 34, NULL, &global_matrix);
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
    MatDestroy(&global_matrix);
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
    MatCreateSeqAIJ(PETSC_COMM_SELF, n_global_dofs, n_global_dofs, 34, NULL, &global_matrix);
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

class Stopwatch {
public:
    inline void start() {
        m_begin = std::chrono::high_resolution_clock::now();
    }
    inline void stop() {
        m_end = std::chrono::high_resolution_clock::now();
    }
    uint microseconds() {
        return std::chrono::duration_cast<std::chrono::microseconds>(m_end - m_begin).count();
    }
    uint milliseconds() {
        return std::chrono::duration_cast<std::chrono::milliseconds>(m_end - m_begin).count();
    }
    uint nanoseconds() {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(m_end - m_begin).count();
    }
private:
    std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::duration<long int, std::ratio<1, 1000000000>>> m_begin;
    std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::duration<long int, std::ratio<1, 1000000000>>> m_end;
};

class TimeWriter {
public:
    TimeWriter(const std::string& filePath) : firstTest(true), file(filePath) {
        file.put('{');
    }
    ~TimeWriter() {
        file.put('}');
        file.put('\n');
    }
    void writeTime(const std::string& testName, double time) {
        if (!firstTest) {
            file.put(',');
        }
        file << '"' << testName << '"' << ": " << time;
        file.flush();
        firstTest = false;
    }
private:
    bool firstTest;
    std::ofstream file;
};

TEST(Spacefilling, space_filling) {

    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    
//     const std::array<std::string, 2> shapes = {"square", "lshape"};
    const std::array<std::string, 1> shapes = {"square"/*, "lshape"*/};
//     const std::array<std::string, 2> structures = {"uniform", "refined"};
    const std::array<std::string, 1> structures = {"uniform"/*, "refined"*/};
//     const std::array<std::string, 2> sizes = {"small", "big"};
    const std::array<std::string, 1> sizes = {/*"small",*/ "big"};
    
    uint progressIndex = 1;
    
    Stopwatch stopwatch;
    
    std::time_t unixTime = std::time(nullptr);
    std::stringstream tmpStream;
    tmpStream << unixTime;
    std::string unixTimeString = tmpStream.str();
    TimeWriter tw("../../../results/multitest_" + unixTimeString + ".json");
    
    for (const std::string& size : sizes) {
        for (const std::string& structure : structures) {
            for (const std::string& shape : shapes) {
                
                const std::string meshName = shape + '_' + structure + '_' + "3D" + '_' + size;
                //optimize_mesh=false is necessary for explicit call of optimizer
                const std::string meshNameJson = " {mesh_file=\"mesh/" + meshName + ".msh\", optimize_mesh=false }";
                
                std::cout << '(' << progressIndex <<  "/40)" << '\n';
                ++progressIndex;
                
                Mesh * mesh = mesh_full_constructor(meshNameJson);
                
                MeshOptimizer backuper(mesh);
                std::vector<Element> elementsBackup = backuper.get_elements(); // method removed
                Armor::Array<double> nodesBackup = backuper.get_nodes(); // method removed
                
                std::cout << meshName << '\n';
                std::cout << "nNodes: " << mesh->n_nodes() << '\n';
                std::cout << "nElements: " << mesh->n_elements() << '\n';
                
                std::array<double, 5> results;
                
//                 {
//                     MeshOptimizer mo(mesh);
//                     std::cout << '(' << progressIndex <<  "/40)" << '\n';
//                     ++progressIndex;
//                     std::cout << "calculating sizes" << '\n';
//                     mo.calculate_sizes();
//                     std::cout << "calculating node curve values" << '\n';
//                     mo.calculate_node_curve_values_as_first_coord();
//                     std::cout << "calculating element curve values" << '\n';
//                     mo.calculate_element_curve_values_as_first_coord();
//                     std::cout << "sorting nodes" << '\n';
//                     mo.sort_nodes();
//                     std::cout << "sorting elements" << '\n';
//                     mo.sort_elements();
//                     std::cout << "getting global matrix" << '\n';
//                     Mat global_matrix = getGlobalMatrix(mesh);
//                     std::cout << "performing calculation" << '\n';
//                     stopwatch.start();
//                     results[1] = performMultiplication(global_matrix, mesh->n_nodes());
//                     stopwatch.stop();
//                     std::cout << "result 2: " << results[1] << '\n';
//                     tw.writeTime(meshName + "_first", stopwatch.microseconds());
//                     std::cout << "copying mesh from backup" << '\n';
//                     mo.set_elements(elementsBackup);
//                     mo.set_nodes(nodesBackup);
//                     MatDestroy(&global_matrix);
//                 }
//                 
//                 {
//                     MeshOptimizer mo(mesh);
//                     std::cout << '(' << progressIndex <<  "/40)" << '\n';
//                     ++progressIndex;
//                     std::cout << "calculating sizes" << '\n';
//                     mo.calculate_sizes();
//                     std::cout << "calculating node curve values" << '\n';
//                     mo.calculate_node_curve_values_as_mean_of_coords();
//                     std::cout << "calculating element curve values" << '\n';
//                     mo.calculate_element_curve_values_as_mean_of_coords();
//                     std::cout << "sorting nodes" << '\n';
//                     mo.sort_nodes();
//                     std::cout << "sorting elements" << '\n';
//                     mo.sort_elements();
//                     std::cout << "getting global matrix" << '\n';
//                     Mat global_matrix = getGlobalMatrix(mesh);
//                     std::cout << "performing calculation" << '\n';
//                     stopwatch.start();
//                     results[2] = performMultiplication(global_matrix, mesh->n_nodes());
//                     stopwatch.stop();
//                     std::cout << "result 3: " << results[2] << '\n';
//                     tw.writeTime(meshName + "_mean", stopwatch.microseconds());
//                     std::cout << "copying mesh from backup" << '\n';
//                     mo.set_elements(elementsBackup);
//                     mo.set_nodes(nodesBackup);
//                     MatDestroy(&global_matrix);
//                 }
//                 
//                 {
//                     MeshOptimizer mo(mesh);
//                     std::cout << '(' << progressIndex <<  "/40)" << '\n';
//                     ++progressIndex;
//                     std::cout << "calculating sizes" << '\n';
//                     mo.calculate_sizes();
//                     std::cout << "calculating node curve values" << '\n';
//                     mo.calculate_node_curve_values_as_zcurve();
//                     std::cout << "calculating element curve values" << '\n';
//                     mo.calculate_element_curve_values_as_zcurve_of_center();
//                     std::cout << "sorting nodes" << '\n';
//                     mo.sort_nodes();
//                     std::cout << "sorting elements" << '\n';
//                     mo.sort_elements();
//                     std::cout << "getting global matrix" << '\n';
//                     Mat global_matrix = getGlobalMatrix(mesh);
//                     std::cout << "performing calculation" << '\n';
//                     stopwatch.start();
//                     results[3] = performMultiplication(global_matrix, mesh->n_nodes());
//                     stopwatch.stop();
//                     std::cout << "result 4: " << results[3] << '\n';
//                     tw.writeTime(meshName + "_zcurve", stopwatch.microseconds());
//                     std::cout << "copying mesh from backup" << '\n';
//                     mo.set_elements(elementsBackup);
//                     mo.set_nodes(nodesBackup);
//                     MatDestroy(&global_matrix);
//                 }
                
                {
                    MeshOptimizer mo(mesh);
                    std::cout << '(' << progressIndex <<  "/40)" << '\n';
                    ++progressIndex;
                    std::cout << "calculating sizes" << '\n';
                    mo.calculate_sizes();
                    std::cout << "calculating node curve values" << '\n';
                    mo.calculate_node_curve_values_as_hilbert();
                    std::cout << "calculating element curve values" << '\n';
                    mo.calculate_element_curve_values_as_hilbert_of_centers();
                    std::cout << "sorting nodes" << '\n';
                    mo.sort_nodes();
                    std::cout << "sorting elements" << '\n';
                    mo.sort_elements();
                    std::cout << "performing calculation" << '\n';
                    stopwatch.start();
                    results[4] = calculation3DAfterSort(mesh);
                    stopwatch.stop();
                    std::cout << "result 5: " << results[4] << '\n';
                    tw.writeTime(meshName + "_hilbert", stopwatch.microseconds());
                    std::cout << "copying mesh from backup" << '\n';
                    mo.set_elements(elementsBackup); // method removed
                    mo.set_nodes(nodesBackup); // method removed
                }
                
                {
                    //
                    //
                    //
                    //
                    //
                    //
                    //
                    //
                    std::cout << "performing calculation" << '\n';
                    stopwatch.start();
                    results[0] = calculation3DBeforeSort(mesh);
                    stopwatch.stop();
                    std::cout << "result 1: " << results[0] << '\n';
                    tw.writeTime(meshName + "_nosort", stopwatch.microseconds());
                }
                
                delete mesh;
                
            }
        }
    }
    
}

