
#ifndef MESH_OPTIMIZER_HH_
#define MESH_OPTIMIZER_HH_

#include <functional>
#include <vector>

#include <armadillo>

#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"
#include "bounding_box.hh"

struct Vertex {
    inline Vertex(arma::vec3 other) {
        for (uint i = 0; i < 3; ++i) {
            data[i] = other[i];
        }
    }
    inline arma::vec3 arma() const {
        return data;
    }
    inline double& operator[](uint i) {
        return data[i];
    }
    inline double operator[](uint i) const {
        return data[i];
    }
    double data[3];
};

struct NodeRef {
    NodeRef(const Vertex& _ref, uint _originalIndex, double _curveValue) : ref(_ref), originalIndex(_originalIndex), curveValue(_curveValue) {}
    std::reference_wrapper<const Vertex> ref;
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

template <uint DIM>
class MeshOptimizer {
    static_assert(DIM == 2 || DIM == 3, "DIM must be either 2 or 3.");
public:
    inline MeshOptimizer(Mesh& _mesh) : mesh(_mesh) {}
    inline void readNodesFromMesh() {
        nodesBackup.reserve(mesh.n_nodes());
        for (uint i = 0; i < mesh.n_nodes(); ++i) {
            nodesBackup.push_back(mesh.nodes_.vec<3>(i));
        }
    }
    inline void readElementsFromMesh() {
        elementsBackup = mesh.element_vec_;
    }
    inline void calculateSizes() {
        nodeSizes.resize(mesh.n_nodes(), INFINITY);
        elementSizes.reserve(mesh.n_elements());
        for (const ElementAccessor<3>& elm : mesh.elements_range()) {
            double elmSize = calculateSizeOfElement(elm);
            elementSizes.push_back(elmSize);
            const Element& el = *elm.element();
            for (uint i = 0; i < DIM + 1; ++i) {
                nodeSizes[el.nodes_[i]] = std::min({nodeSizes[el.nodes_[i]], elmSize});
            }
        }
        std::vector<arma::vec3> tmpArmas;
        tmpArmas.reserve(mesh.n_nodes());
        for (uint i = 0; i < mesh.n_nodes(); ++i) {
            tmpArmas.push_back(nodesBackup[i].arma());
        }
        boundingBox = BoundingBox(tmpArmas);
    }
    inline void calculateNodeCurveValuesAsHilbert() {
        nodeRefs.reserve(mesh.n_nodes());
        for (uint i = 0; i < mesh.n_nodes(); ++i) {
            nodeRefs.emplace_back(nodesBackup[i], i, hilbertValue(normalize(nodesBackup[i], boundingBox), nodeSizes[i]));
        }
    }
    inline void calculateNodeCurveValuesAsZCurve() {
        nodeRefs.reserve(mesh.n_nodes());
        for (uint i = 0; i < mesh.n_nodes(); ++i) {
            nodeRefs.emplace_back(nodesBackup[i], i, zCurveValue(normalize(nodesBackup[i], boundingBox), nodeSizes[i]));
        }
    }
    inline void calculateNodeCurveValuesAsMeanOfCoords() {
        nodeRefs.reserve(mesh.n_nodes());
        for (uint i = 0; i < mesh.n_nodes(); ++i) {
            const Vertex tmpNorm = normalize(nodesBackup[i], boundingBox);
            nodeRefs.emplace_back(nodesBackup[i], i, (tmpNorm[0] + tmpNorm[1] + tmpNorm[2]) / 3);
        }
    }
    inline void calculateNodeCurveValuesAsFirstCoord() {
        nodeRefs.reserve(mesh.n_nodes());
        for (uint i = 0; i < mesh.n_nodes(); ++i) {
            const Vertex tmpNorm = normalize(nodesBackup[i], boundingBox);
            nodeRefs.emplace_back(nodesBackup[i], i, tmpNorm[0]);
        }
    }
    inline void calculateNodeCurveValuesAsObtainedFromElements() {
        nodeRefs.reserve(mesh.n_nodes());
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            const ElementRef& elRf = elementRefs[i];
            for (uint j = 0; j < DIM + 1; ++j) {
                uint nodeIdx = elRf.ref.get().nodes_[j];
                nodeRefs.emplace_back(nodesBackup[nodeIdx], nodeIdx, elementRefs[i].curveValue);
            }
        }
    }
    inline void calculateElementCurveValuesAsMeanOfNodes() {
        elementRefs.reserve(mesh.n_elements());
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            const std::array<uint, 4>& nodeIndexes = elementsBackup[i].nodes_;
            double tmpSum = 0;
            for (uint j = 0; j < DIM + 1; ++j) {
                tmpSum += nodeRefs[nodeIndexes[j]].curveValue;
            }
            elementRefs.emplace_back(elementsBackup[i], tmpSum / DIM + 1);
        }
    }
    inline void calculateElementCurveValuesAsHilbertOfCenters() {
        elementRefs.reserve(mesh.n_elements());
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            elementRefs.emplace_back(elementsBackup[i], hilbertValue(normalize(ElementAccessor<3>(&mesh, i).centre(), boundingBox), elementSizes[i]));
        }
    }
    inline void calculateElementCurveValuesAsZCurveOfCenters() {
        elementRefs.reserve(mesh.n_elements());
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            elementRefs.emplace_back(elementsBackup[i], zCurveValue(normalize(ElementAccessor<3>(&mesh, i).centre(), boundingBox), elementSizes[i]));
        }
    }
    inline void calculateElementCurveValuesAsMeanOfCoords() {
        elementRefs.reserve(mesh.n_elements());
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            arma::vec3 tmpNorm = normalize(ElementAccessor<3>(&mesh, i).centre(), boundingBox);
            nodeRefs.emplace_back(nodesBackup[i], i, (tmpNorm[0] + tmpNorm[1] + tmpNorm[2]) / 3);
        }
    }
    inline void calculateElementCurveValuesAsFirstCoord() {
        elementRefs.reserve(mesh.n_elements());
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            arma::vec3 tmpNorm = normalize(ElementAccessor<3>(&mesh, i).centre(), boundingBox);
            nodeRefs.emplace_back(nodesBackup[i], i, tmpNorm[0]);
        }
    }
    inline void sortNodes() {
        std::sort(nodeRefs.begin(), nodeRefs.end());
        std::vector<uint> newNodeIndexes(nodeRefs.size());
        for (uint i = 0; i < nodeRefs.size(); ++i) {
            newNodeIndexes[nodeRefs[i].originalIndex] = i;
        }
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            for (uint j = 0; j < DIM + 1; ++j) {
                elementsBackup[i].nodes_[j] = newNodeIndexes[elementsBackup[i].nodes_[j]];
            }
        }
    }
    inline void sortElements() {
        std::sort(elementRefs.begin(), elementRefs.end());
    }
    inline void writeNodesToMesh() {
        for (uint i = 0; i < nodeRefs.size(); ++i) {
            mesh.nodes_.set(i) = nodeRefs[i].ref.get().arma();
        }
    }
    inline void writeElementsToMesh() {
        for (uint i = 0; i < elementRefs.size(); ++i) {
            mesh.element_vec_[i] = elementRefs[i].ref.get();
        }
    }
private:
    Mesh& mesh;
    std::vector<Vertex> nodesBackup;
    std::vector<Element> elementsBackup;
    std::vector<NodeRef> nodeRefs;
    std::vector<ElementRef> elementRefs;
    std::vector<double> nodeSizes;
    std::vector<double> elementSizes;
    BoundingBox boundingBox;
    inline double calculateSizeOfElement(const ElementAccessor<3>& elm) {
        return std::min({arma::norm(*elm.node(0) - *elm.node(1)),
                         arma::norm(*elm.node(1) - *elm.node(2)),
                         arma::norm(*elm.node(2) - *elm.node(0)),
                         arma::norm(*elm.node(3) - *elm.node(0)),
                         arma::norm(*elm.node(3) - *elm.node(1)),
                         arma::norm(*elm.node(3) - *elm.node(2))});
    }
    inline double hilbertValue(double x, double y, double eps) {
        if (eps > 1) {
            return 0;
        } else {
            if (x < 0.5) {
                if (y < 0.5) {
                    return hilbertValue(2 * y, 2 * x, 4 * eps) / 4;
                } else {
                    return (1 + hilbertValue(2 * x, 2 * y - 1, 4 * eps)) / 4;
                }
            } else {
                if (y >= 0.5) {
                    return (2 + hilbertValue(2 * x - 1, 2 * y - 1, 4 * eps)) / 4;
                } else {
                    return (3 + hilbertValue(1 - 2 * y, 2 - 2 * x, 4 * eps)) / 4;
                }
            }
        }
    }
    inline double hilbertValue(double x, double y, double z, double eps) {
        if (eps > 1) {
            return 0;
        } else {
            if (z < 0.5) {
                if (x < 0.5) {
                    if (y < 0.5) {
                        return hilbertValue(2 * z, 2 * x, 2 * y, 8 * eps) / 8;
                    } else {
                        return (1 + hilbertValue(2 * y - 1, 2 * z, 2 * x, 8 * eps)) / 8;
                    }
                } else {
                    if (y >= 0.5) {
                        return (2 + hilbertValue(2 * y - 1, 2 * z, 2 * x - 1, 8 * eps)) / 8;
                    } else {
                        return (3 + hilbertValue(2 - 2 * x, 1 - 2 * y, 2 * z, 8 * eps)) / 8;
                    }
                }
            } else {
                if (x >= 0.5) {
                    if (y < 0.5) {
                        return (4 + hilbertValue(2 - 2 * x, 1 - 2 * y, 2 * z - 1, 8 * eps)) / 8;
                    } else {
                        return (5 + hilbertValue(2 * y - 1, 2 - 2 * z, 2 - 2 * x, 8 * eps)) / 8;
                    }
                } else {
                    if (y >= 0.5) {
                        return (6 + hilbertValue(2 * y - 1, 2 - 2 * z, 1 - 2 * x, 8 * eps)) / 8;
                    } else {
                        return (7 + hilbertValue(2 - 2 * z, 2 * x, 1 - 2 * y, 8 * eps)) / 8;
                    }
                }
            }
        }
    }
    inline double zCurveValue(double x, double y, double eps) {
        if (eps > 1) {
            return 0;
        } else {
            if (y < 0.5) {
                if (x < 0.5) {
                    return zCurveValue(2 * x, 2 * y, 4 * eps) / 4;
                } else {
                    return (1 + zCurveValue(2 * x - 1, 2 * y, 4 * eps)) / 4;
                }
            } else {
                if (x < 0.5) {
                    return (2 + zCurveValue(2 * x, 2 * y - 1, 4 * eps)) / 4;
                } else {
                    return (3 + zCurveValue(2 * x - 1, 2 * y - 1, 4 * eps)) / 4;
                }
            }
        }
    }
    inline double zCurveValue(double x, double y, double z, double eps) {
        if (eps > 1) {
            return 0;
        } else {
            if (z < 0.5) {
                if (y < 0.5) {
                    if (x < 0.5) {
                        return zCurveValue(2 * x, 2 * y, 2 * z, 8 * eps) / 8;
                    } else {
                        return (1 + zCurveValue(2 * x - 1, 2 * y, 2 * z, 8 * eps)) / 8;
                    }
                } else {
                    if (x < 0.5) {
                        return (2 + zCurveValue(2 * x, 2 * y - 1, 2 * z, 8 * eps)) / 8;
                    } else {
                        return (3 + zCurveValue(2 * x - 1, 2 * y - 1, 2 * z, 8 * eps)) / 8;
                    }
                }
            } else {
                if (y < 0.5) {
                    if (x < 0.5) {
                        return (4 + zCurveValue(2 * x, 2 * y, 2 * z - 1, 8 * eps)) / 8;
                    } else {
                        return (5 + zCurveValue(2 * x - 1, 2 * y, 2 * z - 1, 8 * eps)) / 8;
                    }
                } else {
                    if (x < 0.5) {
                        return (6 + zCurveValue(2 * x, 2 * y - 1, 2 * z - 1, 8 * eps)) / 8;
                    } else {
                        return (7 + zCurveValue(2 * x - 1, 2 * y - 1, 2 * z - 1, 8 * eps)) / 8;
                    }
                }
            }
        }
    }
    inline double hilbertValue(const Vertex vec, double size) {
        return hilbertValue(vec[0], vec[1], vec[2], size * size * size);
    }
    inline double zCurveValue(const Vertex vec, double size) {
        return zCurveValue(vec[0], vec[1], vec[2], size * size * size);
    }
    inline Vertex normalize(const Vertex& vec, const BoundingBox& boundingBox) {
        const Vertex tmp = arma::vec3(boundingBox.max() - boundingBox.min());
        double scalar = std::max({tmp[0], tmp[1], tmp[2]});
        return arma::vec3((vec.arma() - boundingBox.min()) / scalar);
    }
};

template<>
inline double MeshOptimizer<2>::calculateSizeOfElement(const ElementAccessor<3>& elm) {
    return std::min({arma::norm(*elm.node(0) - *elm.node(1)),
                     arma::norm(*elm.node(1) - *elm.node(2)),
                     arma::norm(*elm.node(2) - *elm.node(0))});
}

template <>
inline double MeshOptimizer<2>::hilbertValue(const Vertex vec, double size) {
    return hilbertValue(vec[0], vec[1], size * size);
}

template <>
inline double MeshOptimizer<2>::zCurveValue(const Vertex vec, double size) {
    return zCurveValue(vec[0], vec[1], size * size);
}

#endif /* MESH_OPTIMIZER_HH_ */
