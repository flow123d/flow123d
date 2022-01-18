/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    output_element.hh
 * @brief   Class OutputElement and its iterator OutputElementIterator on the output mesh.
 */

#ifndef OUTPUT_ELEMENT_HH_
#define OUTPUT_ELEMENT_HH_

#include "element_data_cache.hh"
#include "output_mesh.hh"

template <int spacedim>
class ElementAccessor;


#include "mesh/accessors.hh"
#include "mesh/mesh.h"
#include "mesh/point.hh"
#include "output_mesh.hh"

/** @brief Represents an element of the output mesh.
 * Provides element access on the data of the output mesh (nodes, connectivity, offsets etc.).
 * 
 * Vertex function suppose spacedim = 3 at the moment, hard coded.
 */
class OutputElement
{
public:
    /// Element space dimension = 3.
    static const unsigned int spacedim = 3;
    
    typedef Space<spacedim>::Point Point;
    
    /// Default constructor.
    OutputElement();

    /// Constructor.
    OutputElement(unsigned int ele_idx, std::shared_ptr<OutputMeshBase> output_mesh);
    
    /// @name Output Mesh Getters.
    //@{ 
    unsigned int idx() const;                       ///< Returns index of the output element.
    unsigned int n_nodes() const;                   ///< Returns number of nodes.
    unsigned int dim() const;                       ///< Returns dim of the output element.
    
    /// Returns global index of the node.
    unsigned int node_index(unsigned int loc_idx) const;
    /// Returns global indices of the nodes.
    std::vector<unsigned int> node_list() const;
    Point vertex(unsigned int loc_idx) const;  ///< Returns coordinates of node @p loc_idx.
    std::vector<Point> vertex_list() const;    ///< Returns vector of nodes coordinates.
    
    Point centre() const;                      ///< Computes the barycenter.
    //@}
    
    /// @name Original Mesh Getters.
    //@{ 
    /// Returns pointer to the computational mesh.
    Mesh* orig_mesh() const;
    /// Returns index of the master element in the computational mesh.
    unsigned int orig_element_idx() const;
    ///Gets ElementAccessor of this element
    ElementAccessor<spacedim> element_accessor() const;
    //@}
    
    void inc();

    bool operator==(const OutputElement& other);
private:
    
    /// index of the output element
    unsigned int ele_idx_;
    /// pointer to the output mesh
    std::shared_ptr<OutputMeshBase> output_mesh_;
};

// --------------------------------------------------- OutputElement INLINE implementation -------------------

inline OutputElement::OutputElement() {}

inline OutputElement::OutputElement(unsigned int ele_idx, std::shared_ptr<OutputMeshBase> output_mesh)
: ele_idx_(ele_idx), output_mesh_(output_mesh)
{}

inline void OutputElement::inc()
{
    ele_idx_++;
}

inline bool OutputElement::operator==(const OutputElement& other)
{
    return ele_idx_ == other.ele_idx_;
}


inline Mesh* OutputElement::orig_mesh() const
{
    return output_mesh_->orig_mesh_;
}

inline unsigned int OutputElement::orig_element_idx() const
{
    return (*output_mesh_->orig_element_indices_)[ele_idx_];
}

inline ElementAccessor< OutputElement::spacedim > OutputElement::element_accessor() const
{
    return output_mesh_->orig_mesh_->element_accessor(orig_element_idx());
}


inline unsigned int OutputElement::idx() const
{
    return ele_idx_;
}

inline unsigned int OutputElement::n_nodes() const
{
    return (*output_mesh_->offsets_)[ele_idx_+1] - (*output_mesh_->offsets_)[ele_idx_];
}

inline unsigned int OutputElement::dim() const
{
    return n_nodes()-1;
}


inline unsigned int OutputElement::node_index(unsigned int loc_idx) const
{
    unsigned int n = n_nodes();
    ASSERT(loc_idx < n);
    unsigned int con_off = (*output_mesh_->offsets_)[ele_idx_+1];
    return (* output_mesh_->connectivity_)[con_off - n + loc_idx];
}


inline OutputElement::Point OutputElement::vertex(unsigned int loc_idx) const
{
    unsigned int n = n_nodes();
    ASSERT(loc_idx < n);
    unsigned int con_off = (*output_mesh_->offsets_)[ele_idx_+1];
    unsigned int off = spacedim * (* output_mesh_->connectivity_)[con_off - n + loc_idx];
    auto &d = *( output_mesh_->nodes_->get_data().get() );
    Point point({d[off], d[off+1], d[off+2]});
    return point;
}


inline std::vector< OutputElement::Point > OutputElement::vertex_list() const
{
    const unsigned int n = n_nodes();
    std::vector<Point> vertices(n);
    
    unsigned int con_off = (*output_mesh_->offsets_)[ele_idx_+1];
    auto &d = *( output_mesh_->nodes_->get_data().get() );
    for(unsigned int i=0; i<n; i++) {
        unsigned int off = spacedim * (* output_mesh_->connectivity_)[con_off - n + i];
        vertices[i] = {d[off], d[off+1], d[off+2]};
        off += spacedim;
    }
    return vertices;
}


inline std::vector< unsigned int > OutputElement::node_list() const
{
    unsigned int n = n_nodes();
    unsigned int con_off = (*output_mesh_->offsets_)[ele_idx_+1];
    std::vector<unsigned int> indices(n);
    for(unsigned int i=0; i<n; i++) {
        indices[i] = (* output_mesh_->connectivity_)[con_off - n + i];
    }
    return indices;
}


inline OutputElement::Point OutputElement::centre() const
{
    Point res({0,0,0});
    for(auto& v : vertex_list() ) res += v;
    return res/n_nodes();
}

#endif // OUTPUT_ELEMENT_HH_
