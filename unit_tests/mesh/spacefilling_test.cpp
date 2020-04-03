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
    void setHilbert(std::vector<double>& indices) {
        for (uint i = 0; i < points.size(); ++i) {
            points[i].hilbertIndex = indices[i];
        }
    }
    std::vector<double> getHilbert() {
        std::vector<double> indices(points.size());
        for (uint i = 0; i < points.size(); ++i) {
            indices[i] = points[i].hilbertIndex;
        }
        return indices;
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

double calculationBeforeSort(Mesh * mesh) {
    double checksum = 0;
//     for (uint i = 0; i < 14; ++i) {
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
//     }
    return checksum;
}

double calculationAfterSort(Mesh * mesh) {
    double checksum = 0;
//     for (uint i = 0; i < 14; ++i) {
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
//     }
    return checksum;
}

// class myVec {
// public:
//     myVec() = default;
//     myVec(const arma::vec3& other) {
//         data[0] = other[0];
//         data[1] = other[1];
//         data[2] = other[2];
//     }
//     inline double& operator[](uint i) {
//         return data[i];
//     }
//     inline double operator[](uint i) const {
//         return data[i];
//     }
// private:
//     std::array<double, 3> data;
// };
// 
// bool operator<(const myVec& first, const myVec& second) {
//     if (first[0] != second[0]) {
//         return first[0] < second[0];
//     } else if (first[1] != second[1]) {
//         return first[1] < second[1];
//     } else {
//         return first[2] < second[2];
//     }
// }
// 
// bool operator==(const myVec& first, const myVec& second) {
//     return (first[0] == second[0] && first[1] == second[1] && first[2] == second[2]);
// }

struct nodeRef2D {
    nodeRef2D(const arma::vec3& _ref, uint _originalIndex, double _curveValue) : ref(_ref), originalIndex(_originalIndex), curveValue(_curveValue) {}
    std::reference_wrapper<const arma::vec3> ref;
    uint originalIndex;
    double curveValue;
};

inline bool operator<(const nodeRef2D& first, const nodeRef2D& second) {
    return first.curveValue < second.curveValue;
}

struct elementRef2D {
    elementRef2D(const Element& _ref, double _curveValue) : ref(_ref), curveValue(_curveValue) {}
    std::reference_wrapper<const Element> ref;
    double curveValue;
};

inline bool operator<(const elementRef2D& first, const elementRef2D& second) {
    return first.curveValue < second.curveValue;
}

class MeshOptimalizer2D {
public:
    MeshOptimalizer2D(Mesh& _mesh) : mesh(_mesh) {}
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
        nodeSizes.reserve(mesh.n_nodes());
        elementSizes.reserve(mesh.n_elements());
        for (ElementAccessor<3> elm : mesh.elements_range()) {
            double elmSize = std::min({
                arma::norm(*elm.node(0) - *elm.node(1)),
                arma::norm(*elm.node(1) - *elm.node(2)),
                arma::norm(*elm.node(2) - *elm.node(0))
            });
            elementSizes.push_back(elmSize);
            const Element& el = *elm.element();
            nodeSizes[el.nodes_[0]] = std::min({nodeSizes[el.nodes_[0]], elmSize});
            nodeSizes[el.nodes_[1]] = std::min({nodeSizes[el.nodes_[1]], elmSize});
            nodeSizes[el.nodes_[2]] = std::min({nodeSizes[el.nodes_[2]], elmSize});
        }
    }
    void calculateNodeHilbertValues() {
        nodeRefs.reserve(mesh.n_nodes());
        for (uint i = 0; i < mesh.n_nodes(); ++i) {
            nodeRefs.emplace_back(nodesBackup[i], i, hilbertValue2D(normalize(nodesBackup[i], boundingBox), nodeSizes[i]));
        }
    }
    void calculateElementHibertAsMeanOfNodes() {
        elementRefs.reserve(mesh.n_elements());
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            const std::array<uint, 4>& nodeIndexes = elementRefs[i].ref.get().nodes_;
            elementRefs.emplace_back(elementsBackup[i], 
                                     nodeRefs[nodeIndexes[0]].curveValue 
                                     + nodeRefs[nodeIndexes[1]].curveValue
                                     + nodeRefs[nodeIndexes[2]].curveValue
                                     / 3);
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
        for (uint i = 0; i < elementsBackup.size(); ++i) {
            elementsBackup[i].nodes_[0] = newNodeIndexes[elementsBackup[i].nodes_[0]];
            elementsBackup[i].nodes_[1] = newNodeIndexes[elementsBackup[i].nodes_[1]];
            elementsBackup[i].nodes_[2] = newNodeIndexes[elementsBackup[i].nodes_[2]];
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
            mesh.element_vec_[i] = elementRefs[i].ref;
        }
    }
private:
    Mesh& mesh;
    std::vector<arma::vec3> nodesBackup;
    std::vector<Element> elementsBackup;
    std::vector<nodeRef2D> nodeRefs;
    std::vector<elementRef2D> elementRefs;
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
    double hilbertValue2D(const arma::vec3 vec, double size) {
        return hilbertValue2D(vec[0], vec[1], size * size);
    }
    arma::vec3 normalize(const arma::vec3& vec, const BoundingBox& boundingBox) {
        arma::vec3 tmp = boundingBox.max() - boundingBox.min();
        double scalar = tmp[0] > tmp[1] ? tmp[0] : tmp[1];
        return (vec - boundingBox.min()) / scalar;
    }
};

TEST(Spacefilling, get_centers) {
    Profiler::initialize();

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

//     const std::string testName = "square_uniform;
//     const std::string testName = "square_uniform_2";
//     const std::string testName = "square_refined";
//     const std::string testName = "lshape_refined";
    const std::string testName = "lshape_refined_2";
    const std::string mesh_in_string = "{mesh_file=\"mesh/" + testName + ".msh\"}";
    
    std::cout << mesh_in_string << '\n';
    
    Mesh * mesh = mesh_full_constructor(mesh_in_string);
    
    std::cout << mesh->n_nodes() << '\n';
    
    auto reader = reader_constructor(mesh_in_string);
    reader->read_physical_names(mesh);
    reader->read_raw_mesh(mesh);
    
    std::vector<arma::vec3> centers;
    centers.reserve(mesh->n_elements());
    std::vector<double> elementSizes;
    elementSizes.reserve(mesh->n_elements());
    std::vector<double> nodeSizes(mesh->n_nodes(), INFINITY);
    
    START_TIMER("get_centers_and_calculate_sizes");
    for (ElementAccessor<3> elm : mesh->elements_range()) {
        centers.push_back(elm.centre());
        double elmSize = std::min({
            arma::norm(*elm.node(0) - *elm.node(1)),
            arma::norm(*elm.node(1) - *elm.node(2)),
            arma::norm(*elm.node(2) - *elm.node(0))
        });
        elementSizes.push_back(elmSize);
        const Element& el = *elm.element();
        nodeSizes[el.nodes_[0]] = std::min({nodeSizes[el.nodes_[0]], elmSize});
        nodeSizes[el.nodes_[1]] = std::min({nodeSizes[el.nodes_[1]], elmSize});
        nodeSizes[el.nodes_[2]] = std::min({nodeSizes[el.nodes_[2]], elmSize});
    }
    END_TIMER("get_centers_and_calculate_sizes");

    START_TIMER("calculation_before_sort");
    double checksum1 = calculationBeforeSort(mesh);
    END_TIMER("calculation_before_sort");
    
    std::cout << "checksum 1: " << checksum1 << '\n';
    
    Armor::Array<double>& nodes = mesh->nodes_;
    std::vector<Element>& elements = mesh->element_vec_;
    PointsManager pm;
    
    std::vector<arma::vec3> nodes_backup(mesh->n_nodes());
    for (uint i = 0; i < nodes_backup.size(); ++i) {
        nodes_backup[i] = nodes.vec<3>(i);
    }
    BoundingBox bb(nodes_backup);

    pm.setPoints(nodes_backup);
    pm.getSizes() = nodeSizes;
    pm.setBoundingBox(bb.min()[0], bb.min()[1], bb.max()[0], bb.max()[1]);
    pm.normalize();
    START_TIMER("nodes_hilbert_calculation");
//     pm.calculateHibert(0.0000000000001);
//     pm.calculateHibert(0.0000001);
    pm.calculateHibert();
    END_TIMER("nodes_hilbert_calculation");
    std::vector<double> nodeIndices = pm.getHilbert();
    pm.sortByHilbert();
    std::vector<uint> newNodeIndices = pm.getForwardPermutation();
    
    for (uint i = 0; i < nodes_backup.size(); ++i) {
        nodes.set(newNodeIndices[i]) = nodes_backup[i];
    }
    
    std::vector<double> elementIndices(mesh->n_elements());
    
    for (uint i = 0; i < elementIndices.size(); ++i) {
        elementIndices[i] = (
            nodeIndices[elements[i].nodes_[0]]
            + nodeIndices[elements[i].nodes_[1]]
            + nodeIndices[elements[i].nodes_[2]]
        ) / 3;
    }
    
    vector<Element> elements_backup = elements;
    
    PointsManager pm2;
    pm2.setPoints(centers);
    pm2.getSizes() = elementSizes;
    pm2.setBoundingBox(bb.min()[0], bb.min()[1], bb.max()[0], bb.max()[1]);
    pm2.normalize();
    pm2.calculateHibert();
//     pm2.setHilbert(elementIndices);
    pm2.sortByHilbert();
    std::vector<uint> newElementIndices = pm2.getForwardPermutation();
    for (uint i = 0; i < centers.size(); ++i) {
        elements[newElementIndices[i]] = elements_backup[i];
        elements[i].nodes_[0] = newNodeIndices[elements[i].nodes_[0]];
        elements[i].nodes_[1] = newNodeIndices[elements[i].nodes_[1]];
        elements[i].nodes_[2] = newNodeIndices[elements[i].nodes_[2]];
    }
    
    START_TIMER("calculation_after_sort");
    double checksum2 = calculationAfterSort(mesh);
    END_TIMER("calculation_after_sort");
    
    std::cout << "checksum 2: " << checksum2 << '\n';
    
    // KONTROLA
    
//     std::vector<myVec> nodes_backup_2;
//     nodes_backup_2.reserve(mesh->n_nodes());
//     for (uint i = 0; i < mesh->n_nodes(); ++i) {
//         nodes_backup_2.emplace_back(nodes.vec<3>(i));
//     }
//     
//     std::vector<myVec> nodes_backup_3;
//     nodes_backup_3.reserve(mesh->n_nodes());
//     for (uint i = 0; i < mesh->n_nodes(); ++i) {
//         nodes_backup_3.emplace_back(nodes_backup[i]);
//     }
//     
//     std::sort(nodes_backup_2.begin(), nodes_backup_2.end());
//     std::sort(nodes_backup_3.begin(), nodes_backup_3.end());
//     
//     for (uint i = 0; i < nodes_backup.size(); ++i) {
//         if (!(nodes_backup_2[i] == nodes_backup_3[i])) {
//             std::cout << "wrong i= " << i << '\n';
//             std::cout << "nodes2: " << nodes_backup_2[i][0] << ' ' << nodes_backup_2[i][1] << ' ' << nodes_backup_2[i][2] << '\n';
//             std::cout << "nodes3: " << nodes_backup_3[i][0] << ' ' << nodes_backup_3[i][1] << ' ' << nodes_backup_3[i][2] << '\n' << '\n';
//         }
//     }
//     
//     std::cout << "done" << '\n';
    
//     vector<Element> elements_backup_2 = elements;
//     PointsManager pm3;
//     pm3.setPoints(centers);
//     pm3.getSizes() = elementSizes;
//     pm3.setBoundingBox(bb.min()[0], bb.min()[1], bb.max()[0], bb.max()[1]);
//     pm3.normalize();
//     pm3.calculateHibert();
//     pm3.sortByHilbert();
    
//     std::vector<double> groupIndices;
    
//     pm.fillGroupIndices(groupIndices);
    
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

