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

double calculation2D(Mesh * mesh) {
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

double calculation3D(Mesh * mesh) {
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


class MeshTest : public Mesh {
public:
	enum OptMethod
	{
		hilbert,
		z_curve
	};

	MeshTest(Input::Record in_record, MPI_Comm com = MPI_COMM_WORLD)
    : Mesh(in_record, com) {}

    static MeshTest * mesh_factory(const std::string &input_mesh_str) {
    	Input::Record input_mesh_rec = get_record_accessor(input_mesh_str, Input::FileFormat::format_JSON);
    	MeshTest * mesh = new MeshTest( input_mesh_rec );

    	try {
    		std::shared_ptr< BaseMeshReader > reader = BaseMeshReader::reader_factory(input_mesh_rec.val<FilePath>("mesh_file"));
    		reader->read_physical_names(mesh);
    		reader->read_raw_mesh(mesh);
        } INPUT_CATCH(FilePath::ExcFileOpen, FilePath::EI_Address_String, input_mesh_rec)

        mesh->setup_topology();
        mesh->check_and_finish();
        return mesh;

    }

    Range<ElementAccessor<3>> full_range() const {
    	auto bgn_it = make_iter<ElementAccessor<3>>( ElementAccessor<3>(this, 0) );
    	auto end_it = make_iter<ElementAccessor<3>>( ElementAccessor<3>(this, this->element_vec_.size()) );
    	return Range<ElementAccessor<3>>(bgn_it, end_it);
    }

    inline void optimize_mesh(const std::string &meshName, OptMethod opt_method, Stopwatch &stopwatch)
    {
        std::cout << "optimize mesh" << '\n';
        MeshOptimizer<3> mo(this);
        mo.calculate_sizes();
        if (opt_method == OptMethod::hilbert) {
            mo.calculate_node_curve_values_as_hilbert();
            mo.calculate_element_curve_values_as_hilbert_of_centers();
        } else {
            mo.calculate_node_curve_values_as_zcurve();
            mo.calculate_element_curve_values_as_zcurve_of_center();
        }

        std::cout << "sorting nodes and elements" << '\n';
        this->sort_permuted_nodes_elements( mo.sort_nodes(this->node_permutation_), mo.sort_elements(this->elem_permutation_) );

        std::cout << "performing calculation" << '\n';
        stopwatch.start();
        double result = calculation3D(this);
        stopwatch.stop();
        std::cout << "result: " << result << '\n';
        std::cout << "time: " << stopwatch.microseconds() << '\n';
    }
};


TEST(Spacefilling, space_filling) {

    Profiler::instance();
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    
//     const std::array<std::string, 2> shapes = {"square", "lshape"};
//     const std::array<std::string, 2> structures = {"uniform", "refined"};
//     const std::array<std::string, 2> sizes = {"small", "big"};
    
    uint progressIndex = 1;
    
    Stopwatch stopwatch;
    
//    for (const std::string& size : sizes) {
//        for (const std::string& structure : structures) {
//            for (const std::string& shape : shapes) {
                
                //const std::string meshName = shape + '_' + structure + '_' + "3D" + '_' + size;
                const std::string meshName = "test_27936_elem";
                //optimize_mesh=false is necessary for explicit call of optimizer
                const std::string meshNameJson = " {mesh_file=\"mesh/" + meshName + ".msh\", optimize_mesh=false }";
                
                std::cout << '(' << progressIndex <<  "/40)" << '\n';
                ++progressIndex;
                
                {
                    MeshTest * mesh = MeshTest::mesh_factory(meshNameJson);

                    std::cout << meshName << '\n';
                    std::cout << "nNodes: " << mesh->n_nodes() << '\n';
                    std::cout << "nElements: " << mesh->n_elements() << '\n';

                    mesh->optimize_mesh(meshName, MeshTest::OptMethod::hilbert, stopwatch);
                    delete mesh;
                }

//            }
//        }
//    }
    Profiler::uninitialize();
}

