
#ifndef MESH_OPTIMIZER_HH_
#define MESH_OPTIMIZER_HH_

#include <vector>

#include <armadillo>

#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"
#include "bounding_box.hh"

struct Vec3 {
    inline Vec3(arma::vec3 other) {
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

struct Permutee {
    inline Permutee(uint _originalIndex, double _curveValue) : originalIndex(_originalIndex), curveValue(_curveValue) {}
    uint originalIndex;
    double curveValue;
};

struct Normalizer {
    inline Normalizer() : shift({0, 0, 0}), scalar(1) {}
    inline Normalizer(Vec3 _shift, double _scalar) : shift(_shift), scalar(_scalar) {}
    inline Vec3 normalize(const Vec3 vec) {
        return arma::vec3((vec.arma() - shift.arma()) / scalar);
    }
    inline double normalize(double size) {
        return size / scalar;
    }
    Vec3 shift;
    double scalar;
};

inline bool operator<(const Permutee& first, const Permutee& second) {
    return first.curveValue < second.curveValue;
}

template <uint DIM>
class MeshOptimizer {
    static_assert(DIM == 2 || DIM == 3, "DIM must be either 2 or 3.");
public:
    inline MeshOptimizer(Mesh& _mesh) : mesh(_mesh) {}
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
            tmpArmas.push_back(mesh.nodes_.vec<3>(i));
        }
        BoundingBox boundingBox(tmpArmas);
        const Vec3 dimensions = arma::vec3(boundingBox.max() - boundingBox.min());
        normalizer = Normalizer(boundingBox.min(), std::max({dimensions[0], dimensions[1], dimensions[2]}));
    }
    inline void calculateNodeCurveValuesAsHilbert() {
        nodeRefs.reserve(mesh.n_nodes());
        for (uint i = 0; i < mesh.n_nodes(); ++i) {
            nodeRefs.emplace_back(i, hilbertValue(normalizer.normalize(mesh.nodes_.vec<3>(i)), normalizer.normalize(nodeSizes[i])));
        }
    }
    inline void calculateNodeCurveValuesAsZCurve() {
        nodeRefs.reserve(mesh.n_nodes());
        for (uint i = 0; i < mesh.n_nodes(); ++i) {
            nodeRefs.emplace_back(i, zCurveValue(normalizer.normalize(mesh.nodes_.vec<3>(i)), normalizer.normalize(nodeSizes[i])));
        }
    }
    inline void calculateNodeCurveValuesAsMeanOfCoords() {
        nodeRefs.reserve(mesh.n_nodes());
        for (uint i = 0; i < mesh.n_nodes(); ++i) {
            const Vec3 tmpNorm = normalizer.normalize(mesh.nodes_.vec<3>(i));
            nodeRefs.emplace_back(i, (tmpNorm[0] + tmpNorm[1] + tmpNorm[2]) / 3);
        }
    }
    inline void calculateNodeCurveValuesAsFirstCoord() {
        nodeRefs.reserve(mesh.n_nodes());
        for (uint i = 0; i < mesh.n_nodes(); ++i) {
            const Vec3 tmpNorm = normalizer.normalize(mesh.nodes_.vec<3>(i));
            nodeRefs.emplace_back(i, tmpNorm[0]);
        }
    }
    inline void calculateNodeCurveValuesAsObtainedFromElements() {
        nodeRefs.reserve(mesh.n_nodes());
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            const Element& el = mesh.element_vec_[elementRefs[i].originalIndex];
            for (uint j = 0; j < DIM + 1; ++j) {
                nodeRefs.emplace_back(el.nodes_[j], elementRefs[i].curveValue);
            }
        }
    }
    inline void calculateElementCurveValuesAsMeanOfNodes() {
        elementRefs.reserve(mesh.n_elements());
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            const std::array<uint, 4>& nodeIndexes = mesh.element_vec_[i].nodes_;
            double tmpSum = 0;
            for (uint j = 0; j < DIM + 1; ++j) {
                tmpSum += nodeRefs[nodeIndexes[j]].curveValue;
            }
            elementRefs.emplace_back(i, tmpSum / DIM + 1);
        }
    }
    inline void calculateElementCurveValuesAsHilbertOfCenters() {
        elementRefs.reserve(mesh.n_elements());
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            elementRefs.emplace_back(i, hilbertValue(normalizer.normalize(ElementAccessor<3>(&mesh, i).centre()), normalizer.normalize(elementSizes[i])));
        }
    }
    inline void calculateElementCurveValuesAsZCurveOfCenter() {
        elementRefs.reserve(mesh.n_elements());
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            elementRefs.emplace_back(i, zCurveValue(normalizer.normalize(ElementAccessor<3>(&mesh, i).centre()), normalizer.normalize(elementSizes[i])));
        }
    }
    inline void calculateElementCurveValuesAsMeanOfCoords() {
        elementRefs.reserve(mesh.n_elements());
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            Vec3 tmpNorm = normalizer.normalize(ElementAccessor<3>(&mesh, i).centre());
            elementRefs.emplace_back(i, (tmpNorm[0] + tmpNorm[1] + tmpNorm[2]) / 3);
        }
    }
    inline void calculateElementCurveValuesAsFirstCoord() {
        elementRefs.reserve(mesh.n_elements());
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            Vec3 tmpNorm = normalizer.normalize(ElementAccessor<3>(&mesh, i).centre());
            elementRefs.emplace_back(i, tmpNorm[0]);
        }
    }
    inline void sortNodes() {
        std::sort(nodeRefs.begin(), nodeRefs.end());
        std::vector<uint> newNodeIndexes(nodeRefs.size());
        Armor::Array<double> nodesBackup = mesh.nodes_;
        for (uint i = 0; i < nodeRefs.size(); ++i) {
            newNodeIndexes[nodeRefs[i].originalIndex] = i;
        }
        for (uint i = 0; i < mesh.n_elements(); ++i) {
            for (uint j = 0; j < DIM + 1; ++j) {
                mesh.element_vec_[i].nodes_[j] = newNodeIndexes[mesh.element_vec_[i].nodes_[j]];
            }
        }
        for (uint i = 0; i < nodeRefs.size(); ++i) {
            mesh.nodes_.set(i) = nodesBackup.vec<3>(nodeRefs[i].originalIndex);
        }
    }
    inline void sortElements() {
        std::sort(elementRefs.begin(), elementRefs.end());
        std::vector<Element> elementsBackup = mesh.element_vec_;
        for (uint i = 0; i < elementRefs.size(); ++i) {
            mesh.element_vec_[i] = elementsBackup[elementRefs[i].originalIndex];
        }
    }
    void copyMeshFrom(const Mesh& from) {
        for (uint i = 0; i < from.n_elements(); ++i) {
            for (uint j = 0; j < DIM + 1; ++j) {
                mesh.element_vec_[i].nodes_[j] = from.element_vec_[i].nodes_[j];
            }
        }
        mesh.nodes_ = from.nodes_;
    }
    std::vector<Element> getElements() {
        return mesh.element_vec_;
    }
    Armor::Array<double> getNodes() {
        return mesh.nodes_;
    }
    void setElements(const std::vector<Element>& els) {
        mesh.element_vec_ = els;
    }
    void setNodes(const Armor::Array<double>& nds) {
        mesh.nodes_ = nds;
    }
private:
    Mesh& mesh;
    std::vector<Permutee> nodeRefs;
    std::vector<Permutee> elementRefs;
    std::vector<double> nodeSizes;
    std::vector<double> elementSizes;
    Normalizer normalizer;
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
    inline double hilbertValue(const Vec3 vec, double size) {
        return hilbertValue(vec[0], vec[1], vec[2], size * size * size);
    }
    inline double zCurveValue(const Vec3 vec, double size) {
        return zCurveValue(vec[0], vec[1], vec[2], size * size * size);
    }
};

template<>
inline double MeshOptimizer<2>::calculateSizeOfElement(const ElementAccessor<3>& elm) {
    return std::min({arma::norm(*elm.node(0) - *elm.node(1)),
                     arma::norm(*elm.node(1) - *elm.node(2)),
                     arma::norm(*elm.node(2) - *elm.node(0))});
}

template <>
inline double MeshOptimizer<2>::hilbertValue(const Vec3 vec, double size) {
    return hilbertValue(vec[0], vec[1], size * size);
}

template <>
inline double MeshOptimizer<2>::zCurveValue(const Vec3 vec, double size) {
    return zCurveValue(vec[0], vec[1], size * size);
}

#endif /* MESH_OPTIMIZER_HH_ */
