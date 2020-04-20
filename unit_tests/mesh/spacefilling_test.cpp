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



struct NodeRef {
    NodeRef(const arma::vec3& _ref, uint _originalIndex, double _curveValue) : ref(_ref), originalIndex(_originalIndex), curveValue(_curveValue) {}
    std::reference_wrapper<const arma::vec3> ref;
    uint originalIndex;
    double curveValue;
};

inline bool operator<(const NodeRef& first, const NodeRef& second) {
    return first.curveValue < second.curveValue;
}

struct ElementRef {
    ElementRef(const Element& _ref, double _curveValue) : ref(_ref), curveValue(_curveValue) {}
    std::reference_wrapper<const Element> ref;
    double curveValue;
};

inline bool operator<(const ElementRef& first, const ElementRef& second) {
    return first.curveValue < second.curveValue;
}

class MeshOptimizer {
public:
    MeshOptimizer(Mesh& _mesh) : mesh(_mesh) {}
    void readNodes() {
        nodesBackup.reserve(mesh.n_nodes());
        for (uint i = 0; i < mesh.n_nodes(); ++i) {
            nodesBackup.push_back(mesh.nodes_.vec<3>(i));
        }
        boundingBox = BoundingBox(nodesBackup);
    }
    void readElements() {
        elementsBackup = mesh.element_vec_;
    }
    void calculateSizes() {
        nodeSizes.resize(mesh.n_nodes(), INFINITY);
        elementSizes.reserve(mesh.n_elements());
        for (ElementAccessor<3> elm : mesh.elements_range()) {
            double elmSize = std::min({ // TODO adaptive (element.n_nodes())
                arma::norm(*elm.node(0) - *elm.node(1)),
                arma::norm(*elm.node(1) - *elm.node(2)),
                arma::norm(*elm.node(2) - *elm.node(0))
            });
            elementSizes.push_back(elmSize);
            const Element& el = *elm.element();
            for (uint i = 0; i < el.n_nodes(); ++i) {
                nodeSizes[el.nodes_[i]] = std::min({nodeSizes[el.nodes_[i]], elmSize});
            }
        }
    }
    void calculateNodeHilbertValues2D() {
        nodeRefs.reserve(mesh.n_nodes());
        for (uint i = 0; i < mesh.n_nodes(); ++i) {
            nodeRefs.emplace_back(nodesBackup[i], i, hilbertValue2D(normalize(nodesBackup[i], boundingBox), nodeSizes[i]));
        }
    }
    void calculateNodeHilbertValues3D() {
        nodeRefs.reserve(mesh.n_nodes());
        for (uint i = 0; i < mesh.n_nodes(); ++i) {
            nodeRefs.emplace_back(nodesBackup[i], i, hilbertValue3D(normalize(nodesBackup[i], boundingBox), nodeSizes[i]));
        }
    }
    void calculateElementCurveValueAsMeanOfNodes() {
        elementRefs.reserve(mesh.n_elements());
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            const std::array<uint, 4>& nodeIndexes = elementsBackup[i].nodes_;
            double tmpSum = 0;
            for (uint j = 0; j < elementsBackup[i].n_nodes(); ++j) {
                tmpSum += nodeRefs[nodeIndexes[j]].curveValue;
            }
            elementRefs.emplace_back(elementsBackup[i], tmpSum / elementsBackup[i].n_nodes());
        }
    }
    void calculateElementHibertByCenters() {
        elementRefs.reserve(mesh.n_elements());
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            elementRefs.emplace_back(elementsBackup[i], hilbertValue2D(normalize(ElementAccessor<3>(&mesh, i).centre(), boundingBox), elementSizes[i]));
        }
    }
    void sortNodes() {
        std::sort(nodeRefs.begin(), nodeRefs.end());
        std::vector<uint> newNodeIndexes(nodeRefs.size());
        for (uint i = 0; i < nodeRefs.size(); ++i) {
            newNodeIndexes[nodeRefs[i].originalIndex] = i;
        }
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            for (uint j = 0; j < elementsBackup[i].n_nodes(); ++j) {
                elementsBackup[i].nodes_[j] = newNodeIndexes[elementsBackup[i].nodes_[j]];
            }
        }
    }
    void sortElements() {
        std::sort(elementRefs.begin(), elementRefs.end());
    }
    void exportNodes() {
        for (uint i = 0; i < nodeRefs.size(); ++i) {
            mesh.nodes_.set(i) = nodeRefs[i].ref.get();
        }
    }
    void exportElements() {
        for (uint i = 0; i < elementRefs.size(); ++i) {
            mesh.element_vec_[i] = elementRefs[i].ref.get();
        }
    }
private:
    Mesh& mesh;
    std::vector<arma::vec3> nodesBackup;
    std::vector<Element> elementsBackup;
    std::vector<NodeRef> nodeRefs;
    std::vector<ElementRef> elementRefs;
    std::vector<double> nodeSizes;
    std::vector<double> elementSizes;
    BoundingBox boundingBox;
    double hilbertValue2D(double x, double y, double eps) {
        if (eps > 1) {
            return 0;
        } else {
            if (x < 0.5) {
                if (y < 0.5) {
                    return hilbertValue2D(2 * y, 2 * x, 4 * eps) / 4;
                } else {
                    return (1 + hilbertValue2D(2 * x, 2 * y - 1, 4 * eps)) / 4;
                }
            } else {
                if (y >= 0.5) {
                    return (2 + hilbertValue2D(2 * x - 1, 2 * y - 1, 4 * eps)) / 4;
                } else {
                    return (3 + hilbertValue2D(1 - 2 * y, 2 - 2 * x, 4 * eps)) / 4;
                }
            }
        }
    }
    double hilbertValue3D(double x, double y, double z, double eps) {
        if (eps > 1) {
            return 0;
        } else {
            if (z < 0.5) {
                if (x < 0.5) {
                    if (y < 0.5) {
                        return hilbertValue3D(2 * z, 2 * x, 2 * y, 8 * eps) / 8;
                    } else {
                        return (1 + hilbertValue3D(2 * y - 1, 2 * z, 2 * x, 8 * eps)) / 8;
                    }
                } else {
                    if (y >= 0.5) {
                        return (2 + hilbertValue3D(2 * y - 1, 2 * z, 2 * x - 1, 8 * eps)) / 8;
                    } else {
                        return (3 + hilbertValue3D(2 - 2 * x, 1 - 2 * y, 2 * z, 8 * eps)) / 8;
                    }
                }
            } else {
                if (x >= 0.5) {
                    if (y < 0.5) {
                        return (4 + hilbertValue3D(2 - 2 * x, 1 - 2 * y, 2 * z - 1, 8 * eps)) / 8;
                    } else {
                        return (5 + hilbertValue3D(2 * y - 1, 2 - 2 * z, 2 - 2 * x, 8 * eps)) / 8;
                    }
                } else {
                    if (y >= 0.5) {
                        return (6 + hilbertValue3D(2 * y - 1, 2 - 2 * z, 1 - 2 * x, 8 * eps)) / 8;
                    } else {
                        return (7 + hilbertValue3D(2 - 2 * z, 2 * x, 1 - 2 * y, 8 * eps)) / 8;
                    }
                }
            }
        }
    }
    double hilbertValue2D(const arma::vec3 vec, double size) {
        return hilbertValue2D(vec[0], vec[1], size * size);
    }
    double hilbertValue3D(const arma::vec3 vec, double size) {
        return hilbertValue3D(vec[0], vec[1], vec[2], size * size * size);
    }
    arma::vec3 normalize(const arma::vec3& vec, const BoundingBox& boundingBox) {
        arma::vec3 tmp = boundingBox.max() - boundingBox.min();
        double scalar = std::max({tmp[0], tmp[1], tmp[2]});
        return (vec - boundingBox.min()) / scalar;
    }
};

void printMesh(Mesh& mesh) {
    arma::vec3 sum = {0, 0, 0};
    uint j = 0;
    for (ElementAccessor<3> elm : mesh.elements_range()) {
        const Element& el = *elm.element();
        std::cout << "Element " << j << ":\n";
        for(uint i = 0; i < el.n_nodes(); ++i) {
            arma::vec3 tmp = mesh.nodes_.vec<3>(el.nodes_[i]);
            sum += tmp;
            std::cout << "    " << tmp[0] << ' ' << tmp[1] << ' ' << tmp[2] << '\n';
        }
        ++j;
    }
    std::cout << "SUM = " << sum[0] << ' ' << sum[1] << ' ' << sum[2] << '\n';
}

TEST(Spacefilling, space_filling) {
    Profiler::initialize();

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

//     const std::string testName = "square_uniform;
//     const std::string testName = "square_uniform_2";
//     const std::string testName = "square_refined";
//     const std::string testName = "lshape_refined";
    const std::string testName = "lshape_refined_2_cube";
    const std::string mesh_in_string = "{mesh_file=\"mesh/" + testName + ".msh\"}";
    
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
    
    MeshOptimizer mo(*mesh);
    std::cout << "reading nodes" << '\n';
    mo.readNodes();
    std::cout << "reading elemenst" << '\n';
    mo.readElements();
    std::cout << "calculating sizes" << '\n';
    mo.calculateSizes();
    std::cout << "calculating node hilbert values" << '\n';
    mo.calculateNodeHilbertValues3D();
    std::cout << "calculating element hilbert values" << '\n';
    mo.calculateElementCurveValueAsMeanOfNodes();
    std::cout << "sorting nodes" << '\n';
    mo.sortNodes();
    std::cout << "sorting elements" << '\n';
    mo.sortElements();
    std::cout << "exporting nodes" << '\n';
    mo.exportNodes();
    std::cout << "exporting elements" << '\n';
    mo.exportElements();
    
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
    std::ofstream jsonResult("../../../results/" + testName + "/" + unixTimeString + ".json");
    
    Profiler::instance()->output(jsonResult);
    Profiler::uninitialize();
}

