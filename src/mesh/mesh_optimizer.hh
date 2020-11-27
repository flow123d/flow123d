
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
    inline Permutee(uint original_index, double curve_value) : original_index_(original_index), curve_value_(curve_value) {}
    uint original_index_;
    double curve_value_;
};

struct Normalizer {
    inline Normalizer() : shift_({0, 0, 0}), scalar_(1) {}
    inline Normalizer(Vec3 shift, double scalar) : shift_(shift), scalar_(scalar) {}
    inline Vec3 normalize(const Vec3 vec) {
        return arma::vec3((vec.arma() - shift_.arma()) / scalar_);
    }
    inline double normalize(double size) {
        return size / scalar_;
    }
    Vec3 shift_;
    double scalar_;
};

inline bool operator<(const Permutee& first, const Permutee& second) {
    return first.curve_value_ < second.curve_value_;
}

class MeshOptimizer {
public:
    inline MeshOptimizer(Mesh* mesh) : mesh_(mesh) {}

    inline void calculate_sizes() {
        node_sizes_.resize(mesh_->n_nodes(), INFINITY);
        element_sizes_.reserve(mesh_->n_elements());
        for (const ElementAccessor<3>& elm : mesh_->elements_range()) {
            double elm_size = this->calculate_size_of_element(elm);
            element_sizes_.push_back(elm_size);
            const Element& el = *elm.element();
            for (uint i = 0; i < elm.dim() + 1; ++i) {
                node_sizes_[el.nodes_[i]] = std::min({node_sizes_[el.nodes_[i]], elm_size});
            }
        }
        std::vector<arma::vec3> tmp_armas;
        tmp_armas.reserve(mesh_->n_nodes());
        for (uint i = 0; i < mesh_->n_nodes(); ++i) {
            tmp_armas.push_back(mesh_->nodes_.vec<3>(i));
        }
        BoundingBox bounding_box(tmp_armas);
        normalizer_ = Normalizer(bounding_box.min(), arma::vec3(bounding_box.max() - bounding_box.min()).max());
    }

    inline void calculate_node_curve_values_as_hilbert() {
        node_refs_.reserve(mesh_->n_nodes());
        for (uint i = 0; i < mesh_->n_nodes(); ++i) {
            node_refs_.emplace_back(i, hilbert_value(normalizer_.normalize(mesh_->nodes_.vec<3>(i)), normalizer_.normalize(node_sizes_[i]), mesh_->element_vec_[i].dim()));
        }
    }

    inline void calculate_node_curve_values_as_zcurve() {
        node_refs_.reserve(mesh_->n_nodes());
        for (uint i = 0; i < mesh_->n_nodes(); ++i) {
            node_refs_.emplace_back(i, zcurve_value(normalizer_.normalize(mesh_->nodes_.vec<3>(i)), normalizer_.normalize(node_sizes_[i]), mesh_->element_vec_[i].dim()));
        }
    }

    inline void calculate_node_curve_values_as_obtained_from_elements() {
        node_refs_.reserve(mesh_->n_nodes());
        for (uint i = 0; i < mesh_->n_elements(); ++i) {
            const Element& el = mesh_->element_vec_[element_refs_[i].original_index_];
            for (uint j = 0; j < el.dim() + 1; ++j) {
                node_refs_.emplace_back(el.nodes_[j], element_refs_[i].curve_value_);
            }
        }
    }

    inline void calculate_element_curve_values_as_hilbert_of_centers() {
        element_refs_.reserve(mesh_->n_elements());
        for (uint i = 0; i < mesh_->n_elements(); ++i) {
            element_refs_.emplace_back(i, hilbert_value(normalizer_.normalize(ElementAccessor<3>(mesh_, i).centre()), normalizer_.normalize(element_sizes_[i]), mesh_->element_vec_[i].dim()));
        }
    }

    inline void calculate_element_curve_values_as_zcurve_of_center() {
        element_refs_.reserve(mesh_->n_elements());
        for (uint i = 0; i < mesh_->n_elements(); ++i) {
            element_refs_.emplace_back(i, zcurve_value(normalizer_.normalize(ElementAccessor<3>(mesh_, i).centre()), normalizer_.normalize(element_sizes_[i]), mesh_->element_vec_[i].dim()));
        }
    }

    inline std::vector<int> sort_nodes(std::vector<unsigned int> & node_permutation) {
        return this->sort(node_refs_, node_permutation);
    }

    inline std::vector<int> sort_elements(std::vector<unsigned int> & elem_permutation) {
        return this->sort(element_refs_, elem_permutation);
    }

private:
    Mesh* mesh_;
    std::vector<Permutee> node_refs_;
    std::vector<Permutee> element_refs_;
    std::vector<double> node_sizes_;
    std::vector<double> element_sizes_;
    Normalizer normalizer_;

    inline double calculate_size_of_element(const ElementAccessor<3>& elm) const {
    	switch (elm.dim()) {
    	case 1:
    	    return arma::norm(*elm.node(0) - *elm.node(1));
    	case 2:
    	    return std::min({arma::norm(*elm.node(0) - *elm.node(1)),
                             arma::norm(*elm.node(1) - *elm.node(2)),
                             arma::norm(*elm.node(2) - *elm.node(0))});
    	case 3:
    	    return std::min({arma::norm(*elm.node(0) - *elm.node(1)),
                             arma::norm(*elm.node(1) - *elm.node(2)),
                             arma::norm(*elm.node(2) - *elm.node(0)),
                             arma::norm(*elm.node(3) - *elm.node(0)),
                             arma::norm(*elm.node(3) - *elm.node(1)),
                             arma::norm(*elm.node(3) - *elm.node(2))});
    	default:
    	    return 0.0; // should not happen
    	}
    }

    inline double hilbert_value(double x, double eps) const {
        if (eps > 1) {
            return 0;
        } else {
            if (x < 0.5) {
                return hilbert_value(2 * x, 2 * eps) / 2;
            } else {
                return (1 + hilbert_value(2 * x - 1, 2 * eps)) / 2;
            }
        }
    }

    inline double hilbert_value(double x, double y, double eps) const {
        if (eps > 1) {
            return 0;
        } else {
            if (x < 0.5) {
                if (y < 0.5) {
                    return hilbert_value(2 * y, 2 * x, 4 * eps) / 4;
                } else {
                    return (1 + hilbert_value(2 * x, 2 * y - 1, 4 * eps)) / 4;
                }
            } else {
                if (y >= 0.5) {
                    return (2 + hilbert_value(2 * x - 1, 2 * y - 1, 4 * eps)) / 4;
                } else {
                    return (3 + hilbert_value(1 - 2 * y, 2 - 2 * x, 4 * eps)) / 4;
                }
            }
        }
    }

    inline double hilbert_value(double x, double y, double z, double eps) const {
        if (eps > 1) {
            return 0;
        } else {
            if (z < 0.5) {
                if (x < 0.5) {
                    if (y < 0.5) {
                        return hilbert_value(2 * z, 2 * x, 2 * y, 8 * eps) / 8;
                    } else {
                        return (1 + hilbert_value(2 * y - 1, 2 * z, 2 * x, 8 * eps)) / 8;
                    }
                } else {
                    if (y >= 0.5) {
                        return (2 + hilbert_value(2 * y - 1, 2 * z, 2 * x - 1, 8 * eps)) / 8;
                    } else {
                        return (3 + hilbert_value(2 - 2 * x, 1 - 2 * y, 2 * z, 8 * eps)) / 8;
                    }
                }
            } else {
                if (x >= 0.5) {
                    if (y < 0.5) {
                        return (4 + hilbert_value(2 - 2 * x, 1 - 2 * y, 2 * z - 1, 8 * eps)) / 8;
                    } else {
                        return (5 + hilbert_value(2 * y - 1, 2 - 2 * z, 2 - 2 * x, 8 * eps)) / 8;
                    }
                } else {
                    if (y >= 0.5) {
                        return (6 + hilbert_value(2 * y - 1, 2 - 2 * z, 1 - 2 * x, 8 * eps)) / 8;
                    } else {
                        return (7 + hilbert_value(2 - 2 * z, 2 * x, 1 - 2 * y, 8 * eps)) / 8;
                    }
                }
            }
        }
    }

    inline double zcurve_value(double x, double eps) const {
        if (eps > 1) {
            return 0;
        } else {
            if (x < 0.5) {
                return zcurve_value(2 * x, 2 * eps) / 2;
            } else {
                return (1 + zcurve_value(2 * x - 1, 2 * eps)) / 2;
            }
        }
    }

    inline double zcurve_value(double x, double y, double eps) const {
        if (eps > 1) {
            return 0;
        } else {
            if (y < 0.5) {
                if (x < 0.5) {
                    return zcurve_value(2 * x, 2 * y, 4 * eps) / 4;
                } else {
                    return (1 + zcurve_value(2 * x - 1, 2 * y, 4 * eps)) / 4;
                }
            } else {
                if (x < 0.5) {
                    return (2 + zcurve_value(2 * x, 2 * y - 1, 4 * eps)) / 4;
                } else {
                    return (3 + zcurve_value(2 * x - 1, 2 * y - 1, 4 * eps)) / 4;
                }
            }
        }
    }

    inline double zcurve_value(double x, double y, double z, double eps) const {
        if (eps > 1) {
            return 0;
        } else {
            if (z < 0.5) {
                if (y < 0.5) {
                    if (x < 0.5) {
                        return zcurve_value(2 * x, 2 * y, 2 * z, 8 * eps) / 8;
                    } else {
                        return (1 + zcurve_value(2 * x - 1, 2 * y, 2 * z, 8 * eps)) / 8;
                    }
                } else {
                    if (x < 0.5) {
                        return (2 + zcurve_value(2 * x, 2 * y - 1, 2 * z, 8 * eps)) / 8;
                    } else {
                        return (3 + zcurve_value(2 * x - 1, 2 * y - 1, 2 * z, 8 * eps)) / 8;
                    }
                }
            } else {
                if (y < 0.5) {
                    if (x < 0.5) {
                        return (4 + zcurve_value(2 * x, 2 * y, 2 * z - 1, 8 * eps)) / 8;
                    } else {
                        return (5 + zcurve_value(2 * x - 1, 2 * y, 2 * z - 1, 8 * eps)) / 8;
                    }
                } else {
                    if (x < 0.5) {
                        return (6 + zcurve_value(2 * x, 2 * y - 1, 2 * z - 1, 8 * eps)) / 8;
                    } else {
                        return (7 + zcurve_value(2 * x - 1, 2 * y - 1, 2 * z - 1, 8 * eps)) / 8;
                    }
                }
            }
        }
    }

    inline double hilbert_value(const Vec3 vec, double size, unsigned int dim) const {
        switch (dim) {
        case 1:
            return hilbert_value(vec[0], size);
        case 2:
            return hilbert_value(vec[0], vec[1], size * size);
        case 3:
            return hilbert_value(vec[0], vec[1], vec[2], size * size * size);
        default:
            return 0.0; // should not happen
        }
    }

    inline double zcurve_value(const Vec3 vec, double size, unsigned int dim) const {
        switch (dim) {
        case 1:
            return zcurve_value(vec[0], size);
        case 2:
            return zcurve_value(vec[0], vec[1], size * size);
        case 3:
            return zcurve_value(vec[0], vec[1], vec[2], size * size * size);
        default:
            return 0.0; // should not happen
        }
    }

    inline std::vector<int> sort(std::vector<Permutee> &refs, std::vector<unsigned int> &mesh_perm) {
    	ASSERT_DBG(refs.size() <= mesh_perm.size());
        std::sort(refs.begin(), refs.end());
        std::vector<int> mesh_ids;
        mesh_ids.reserve(refs.size());
        for (uint i = 0; i < refs.size(); ++i) {
            mesh_perm[ refs[i].original_index_ ] = i;
            mesh_ids.push_back( refs[i].original_index_ );
        }
        return mesh_ids;
    }

};

/*template<>
inline double MeshOptimizer<2>::calculate_size_of_element(const ElementAccessor<3>& elm) {
    return std::min({arma::norm(*elm.node(0) - *elm.node(1)),
                     arma::norm(*elm.node(1) - *elm.node(2)),
                     arma::norm(*elm.node(2) - *elm.node(0))});
}

template <>
inline double MeshOptimizer<2>::hilbert_value(const Vec3 vec, double size) {
    return hilbert_value(vec[0], vec[1], size * size);
}

template <>
inline double MeshOptimizer<2>::zcurve_value(const Vec3 vec, double size) {
    return zcurve_value(vec[0], vec[1], size * size);
}*/

#endif /* MESH_OPTIMIZER_HH_ */
