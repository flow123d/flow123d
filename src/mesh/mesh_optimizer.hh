
#ifndef MESH_OPTIMIZER_HH_
#define MESH_OPTIMIZER_HH_

#include <vector>

#include <armadillo>

#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"
#include "bounding_box.hh"

struct Permutee {
    inline Permutee(uint original_index, double curve_value) : original_index_(original_index), curve_value_(curve_value) {}
    uint original_index_;
    double curve_value_;
};

struct Normalizer {
    inline Normalizer() : shift_({0, 0, 0}), scalar_(1) {}
    inline Normalizer(arma::vec3 shift, double scalar) : shift_(shift), scalar_(scalar) {}
    inline arma::vec3 normalize(const arma::vec3 vec) {
        return (vec - shift_) / scalar_;
    }
    arma::vec3 shift_;
    double scalar_;
};

inline bool operator<(const Permutee& first, const Permutee& second) {
    return first.curve_value_ < second.curve_value_;
}

template <uint DIM>
class MeshOptimizer {
    static_assert(DIM == 2 || DIM == 3, "DIM must be either 2 or 3.");
public:
    inline MeshOptimizer(Mesh* mesh) : mesh_(mesh) {}

    inline void calculate_sizes() {
    	BoundingBox bounding_box( *mesh_->node(0) );
    	for (auto nod : mesh_->node_range()) {
    	    bounding_box.expand(*nod);
    	}
    	double mesh_longest_size = bounding_box.longest_size();
    	normalizer_ = Normalizer(bounding_box.min(), mesh_longest_size);

        node_sizes_.resize(mesh_->n_nodes(), INFINITY);
        element_sizes_.reserve(mesh_->n_elements());
        for (const ElementAccessor<3>& elm : mesh_->elements_range()) {
            double elm_norm_size = elm.bounding_box().longest_size() / mesh_longest_size;
            element_sizes_.push_back(elm_norm_size);
            for (uint i = 0; i < elm.dim() + 1; ++i) {
                node_sizes_[elm->node_idx(i)] = std::min({ node_sizes_[elm->node_idx(i)], elm_norm_size });
            }
        }
    }

    inline void calculate_node_curve_values_as_hilbert() {
        node_refs_.reserve(mesh_->n_nodes());
        uint i = 0;
        for (auto nod : mesh_->node_range()) {
            node_refs_.emplace_back(i, hilbert_value(normalizer_.normalize( *nod ), node_sizes_[i]));
            i++;
        }
    }

    inline void calculate_node_curve_values_as_zcurve() {
        node_refs_.reserve(mesh_->n_nodes());
        uint i = 0;
        for (auto nod : mesh_->node_range()) {
            node_refs_.emplace_back(i, zcurve_value(normalizer_.normalize( *nod ), node_sizes_[i]));
        }
    }

    inline void calculate_node_curve_values_as_obtained_from_elements() {
        node_refs_.reserve(mesh_->n_nodes());
        for (uint i = 0; i < mesh_->n_elements(); ++i) {
        	auto elm = mesh_->element_accessor(element_refs_[i].original_index_);
            for (uint j = 0; j < elm->dim() + 1; ++j) {
                node_refs_.emplace_back(elm->node_idx(j), element_refs_[i].curve_value_);
            }
        }
    }

    inline void calculate_element_curve_values_as_hilbert_of_centers() {
        element_refs_.reserve(mesh_->n_elements());
        uint i = 0;
        for (auto elm : mesh_->elements_range()) {
            element_refs_.emplace_back(i, hilbert_value(normalizer_.normalize(elm.centre()), element_sizes_[i]));
            i++;
        }
    }

    inline void calculate_element_curve_values_as_zcurve_of_center() {
        element_refs_.reserve(mesh_->n_elements());
        uint i = 0;
        for (auto elm : mesh_->elements_range()) {
            element_refs_.emplace_back(i, zcurve_value(normalizer_.normalize(elm.centre()), element_sizes_[i]));
            i++;
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

    inline double hilbert_value(arma::vec3 vec, double size) const {
        return hilbert_value(vec[0], vec[1], vec[2], size * size * size);
    }

    inline double zcurve_value(arma::vec3 vec, double size) const {
        return zcurve_value(vec[0], vec[1], vec[2], size * size * size);
    }

    inline std::vector<int> sort(std::vector<Permutee> &refs, std::vector<unsigned int> &mesh_perm) {
    	ASSERT(refs.size() <= mesh_perm.size());
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

template<>
inline double MeshOptimizer<2>::hilbert_value(arma::vec3 vec, double size) const {
    return hilbert_value(vec[0], vec[1], size * size);
}

template <>
inline double MeshOptimizer<2>::zcurve_value(arma::vec3 vec, double size) const {
    return zcurve_value(vec[0], vec[1], size * size);
}

#endif /* MESH_OPTIMIZER_HH_ */
