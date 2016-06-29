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
 * Vertex function suppose spacedim = 3 at the moment.
 */
class OutputElement
{
public:
    /// Element space dimension = 3.
    static const unsigned int spacedim = 3;
    
    typedef Space<spacedim>::Point Point;
    
    /// Constructor.
    OutputElement(unsigned int ele_idx, OutputMesh* output_mesh);
    
    /// @name Output Mesh Getters.
    //@{ 
    unsigned int idx() const;                       ///< Returns index of the output element.
    unsigned int n_nodes() const;                   ///< Returns number of nodes.
    unsigned int dim() const;                       ///< Returns dim of the output element.
    
    /// Returns global index of the node.
    unsigned int node_index(unsigned int loc_idx) const;
    Point vertex(unsigned int loc_idx) const;  ///< Returns coordinates of node @p loc_idx.
    std::vector<Point> vertex_list() const;    ///< Returns vector of nodes coordinates.
    
    /// Returns global index of the node. (DISCONTINOUS)
    unsigned int node_index_disc(unsigned int loc_idx) const;
    
    //TODO: vertex_disc and vertex_list_disc return the same coordinates as vertex and vertex_list
    // It is meaningful only if we have ONLY discontinuous vectors for nodes and connectivity
    /// Returns coordinates of node @p loc_idx. (DISCONTINOUS)
    Point vertex_disc(unsigned int loc_idx) const;
    /// Returns vector of nodes coordinates. (DISCONTINOUS)
    std::vector<Point> vertex_list_disc() const;
    
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
    
    void operator++();
    bool operator==(const OutputElement& other);
private:
    
    /// Returns global index of the node. (DISCONTINUOUS)
    unsigned int node_index_internal(unsigned int loc_idx, 
                                     std::shared_ptr<MeshData<unsigned int>> connectivity) const;
    /// Returns coordinates of node @p loc_idx.  (DISCONTINUOUS)
    Point vertex_internal(unsigned int loc_idx,
                               std::shared_ptr<MeshData<unsigned int>> connectivity,
                               std::shared_ptr<MeshData<double>> nodes) const;
    /// Returns vector of nodes coordinates.  (DISCONTINUOUS)
    std::vector<Point> vertex_list_internal(std::shared_ptr<MeshData<unsigned int>> connectivity,
                                            std::shared_ptr<MeshData<double>> nodes) const;
    
//     friend void OutputElementIterator::operator++();
    unsigned int ele_idx_;      ///< index of the output element
    OutputMesh* output_mesh_;   ///< pointer to the output mesh
};

// --------------------------------------------------- OutputElement INLINE implementation -------------------

inline OutputElement::OutputElement(unsigned int ele_idx, OutputMesh* output_mesh)
: ele_idx_(ele_idx), output_mesh_(output_mesh)
{}

inline void OutputElement::operator++()
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
    if(ele_idx_ == 0)
        return (*output_mesh_->offsets_)[0];
    else
        return (*output_mesh_->offsets_)[ele_idx_] - (*output_mesh_->offsets_)[ele_idx_-1];
}

inline unsigned int OutputElement::dim() const
{
    return n_nodes()-1;
}

inline unsigned int OutputElement::node_index_internal(unsigned int loc_idx,
                                                shared_ptr< MeshData< unsigned int > > connectivity) const
{
    unsigned int n = n_nodes();
    ASSERT_DBG(loc_idx < n);
    unsigned int con_off = (*output_mesh_->offsets_)[ele_idx_];
    return (*connectivity)[con_off - n + loc_idx];
}

inline OutputElement::Point OutputElement::vertex_internal(unsigned int loc_idx,
                                                 shared_ptr< MeshData< unsigned int > > connectivity,
                                                 shared_ptr< MeshData< double > > nodes) const
{
    unsigned int n = n_nodes();
    ASSERT_DBG(loc_idx < n);
    unsigned int con_off = (*output_mesh_->offsets_)[ele_idx_];
    unsigned int off = spacedim * (*connectivity)[con_off - n + loc_idx];
    auto &d = nodes->data_;
    Point point({d[off], d[off+1], d[off+2]});
    return point;
}

inline std::vector< OutputElement::Point > OutputElement::vertex_list_internal(shared_ptr< MeshData< unsigned int > > connectivity, 
                                                                               shared_ptr< MeshData< double > > nodes) const
{
    const unsigned int n = n_nodes();
    std::vector<Point> vertices(n);
    
    unsigned int con_off = (*output_mesh_->offsets_)[ele_idx_];
    auto &d = nodes->data_;
    for(unsigned int i=0; i<n; i++) {
        unsigned int off = spacedim * (*connectivity)[con_off - n + i];
        vertices[i] = {d[off], d[off+1], d[off+2]};
        off += spacedim;
    }
    return vertices;
}


inline unsigned int OutputElement::node_index(unsigned int loc_idx) const
{
    return node_index_internal(loc_idx, output_mesh_->connectivity_);
}

inline unsigned int OutputElement::node_index_disc(unsigned int loc_idx) const
{
    ASSERT_PTR(output_mesh_->discont_connectivity_);
    return node_index_internal(loc_idx, output_mesh_->discont_connectivity_);
}

inline OutputElement::Point OutputElement::vertex(unsigned int loc_idx) const
{
    return vertex_internal(loc_idx, output_mesh_->connectivity_, output_mesh_->nodes_);
}

inline OutputElement::Point OutputElement::vertex_disc(unsigned int loc_idx) const
{
    return vertex_internal(loc_idx, output_mesh_->discont_connectivity_, output_mesh_->discont_nodes_);
}

inline std::vector< OutputElement::Point > OutputElement::vertex_list() const
{
    return vertex_list_internal(output_mesh_->connectivity_, output_mesh_->nodes_);
}

inline std::vector< OutputElement::Point > OutputElement::vertex_list_disc() const
{
    return vertex_list_internal(output_mesh_->discont_connectivity_, output_mesh_->discont_nodes_);
}

inline OutputElement::Point OutputElement::centre() const
{
    Point res({0,0,0});
    for(auto& v : vertex_list() ) res += v;
    return res/n_nodes();
}

#endif // OUTPUT_ELEMENT_HH_