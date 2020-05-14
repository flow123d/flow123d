/*!
 *
 * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
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
 * @file    accessors_impl.hh
 * @brief   Implementation of the inline functions of the mesh accessors.
 */

/*******************************************************************************
 * ElementAccessor IMPLEMENTATION
 *******************************************************************************/
/**
 * Default invalid accessor.
 */
template <int spacedim> inline
ElementAccessor<spacedim>::ElementAccessor()
: mesh_(nullptr)
{}

/**
 * Regional accessor.
 */
template <int spacedim> inline
ElementAccessor<spacedim>::ElementAccessor(const Mesh *mesh, RegionIdx r_idx)
: dim_(undefined_dim_),
  mesh_(mesh),
  r_idx_(r_idx)
{}

/**
 * Element accessor.
 */
template <int spacedim> inline
ElementAccessor<spacedim>::ElementAccessor(const Mesh *mesh, unsigned int idx)
: mesh_(mesh),
  boundary_(idx>=mesh->n_elements()),
  element_idx_(idx),
  r_idx_(element()->region_idx())
{
    dim_=element()->dim();
}

template <int spacedim> inline
void ElementAccessor<spacedim>::inc() {
    ASSERT(!is_regional()).error("Do not call inc() for regional accessor!");
    element_idx_++;
    r_idx_ = element()->region_idx();
    dim_=element()->dim();
    boundary_ = (element_idx_>=mesh_->n_elements());
}

template <int spacedim> inline
vector<arma::vec3> ElementAccessor<spacedim>::vertex_list() const {
    vector<arma::vec3> vertices(element()->n_nodes());
    for(unsigned int i=0; i<element()->n_nodes(); i++) vertices[i]=*node(i);
    return vertices;
}

template <int spacedim> inline
double ElementAccessor<spacedim>::tetrahedron_jacobian() const
{
    ASSERT(dim() == 3)(dim()).error("Cannot provide Jacobian for dimension other than 3.");
    return arma::dot( arma::cross(*( node(1) ) - *( node(0) ),
                                    *( node(2) ) - *( node(0) )),
                    *( node(3) ) - *( node(0) )
                    );
}

/**
 * SET THE "METRICS" FIELD IN STRUCT ELEMENT
 */
template <int spacedim> inline
double ElementAccessor<spacedim>::measure() const {
    switch (dim()) {
        case 0:
            return 1.0;
            break;
        case 1:
            return arma::norm(*( node(1) ) - *( node(0) ) , 2);
            break;
        case 2:
            return
                arma::norm(
                    arma::cross(*( node(1) ) - *( node(0) ), *( node(2) ) - *( node(0) )),
                    2
                ) / 2.0 ;
            break;
        case 3:
            return fabs(
                arma::dot(
                    arma::cross(*( node(1) ) - *( node(0) ), *( node(2) ) - *( node(0) )),
                    *( node(3) ) - *( node(0) ) )
                ) / 6.0;
            break;
    }
    return 1.0;
}

/**
 * SET THE "CENTRE[]" FIELD IN STRUCT ELEMENT
 */
template <int spacedim> inline
arma::vec::fixed<spacedim> ElementAccessor<spacedim>::centre() const {
	ASSERT(is_valid()).error("Invalid element accessor.");
    if (is_regional() ) return arma::vec::fixed<spacedim>();

    arma::vec::fixed<spacedim> centre;
    centre.zeros();

    for (unsigned int li=0; li<element()->n_nodes(); li++) {
        centre += *node( li );
    }
    centre /= (double) element()->n_nodes();
    return centre;
}


template <int spacedim> inline
double ElementAccessor<spacedim>::quality_measure_smooth(SideIter side) const {
    if (dim_==3) {
        double sum_faces=0;
        double face[4];
        for(unsigned int i=0; i<4; i++, ++side) sum_faces+=( face[i]=side->measure());

        double sum_pairs=0;
        for(unsigned int i=0;i<3;i++)
            for(unsigned int j=i+1;j<4;j++) {
                unsigned int i_line = RefElement<3>::line_between_faces(i,j);
                arma::vec line = *node(RefElement<3>::interact(Interaction<0,1>(i_line))[1]) - *node(RefElement<3>::interact(Interaction<0,1>(i_line))[0]);
                sum_pairs += face[i]*face[j]*arma::dot(line, line);
            }
        double regular = (2.0*sqrt(2.0/3.0)/9.0); // regular tetrahedron
        return fabs( measure()
                * pow( sum_faces/sum_pairs, 3.0/4.0))/ regular;

    }
    if (dim_==2) {
        return fabs(
                measure()/
                pow(
                         arma::norm(*node(1) - *node(0), 2)
                        *arma::norm(*node(2) - *node(1), 2)
                        *arma::norm(*node(0) - *node(2), 2)
                        , 2.0/3.0)
               ) / ( sqrt(3.0) / 4.0 ); // regular triangle
    }
    return 1.0;
}

template <int spacedim> inline
SideIter ElementAccessor<spacedim>::side(const unsigned int loc_index) {
    return SideIter( Side(mesh_, element_idx_, loc_index) );
}

template <int spacedim> inline
const SideIter ElementAccessor<spacedim>::side(const unsigned int loc_index) const {
    return SideIter( Side(mesh_, element_idx_, loc_index) );
}



/*******************************************************************************
 * Edge IMPLEMENTATION
 *******************************************************************************/

inline Edge::Edge()
: mesh_(nullptr),
  edge_idx_(Mesh::undef_idx)
{}

inline Edge::Edge(const Mesh *mesh, unsigned int edge_idx)
: mesh_(mesh),
  edge_idx_(edge_idx)
{}

inline const EdgeData* Edge::edge_data() const
{
    ASSERT_DBG(is_valid());
    ASSERT_LT_DBG(edge_idx_, mesh_->edges.size());
    return &mesh_->edges[edge_idx_];
}

inline SideIter Edge::side(const unsigned int i) const {
    return edge_data()->side_[i];
}




/*******************************************************************************
 * Side IMPLEMENTATION
 *******************************************************************************/

inline Side::Side()
: mesh_(NULL), elem_idx_(0), side_idx_(0)
{}

inline Side::Side(const Mesh * mesh, unsigned int elem_idx, unsigned int set_lnum)
: mesh_(mesh), elem_idx_(elem_idx), side_idx_(set_lnum)
{
	mesh_->check_element_size(elem_idx);
}

inline unsigned int Side::dim() const {
        return element()->dim()-1;
    }

// returns true for all sides either on boundary or connected to vb neigboring
inline bool Side::is_external() const {
    return edge().n_sides() == 1;
}

// returns true for all sides either on boundary or connected to vb neigboring
inline bool Side::is_boundary() const {
    return is_external() && cond_idx() != Mesh::undef_idx;
}

inline NodeAccessor<3> Side::node(unsigned int i) const {
    int i_n = mesh_->side_nodes[dim()][side_idx_][i];

    return element().node( i_n );
}

inline ElementAccessor<3> Side::element() const {
    ASSERT( is_valid() ).error("Wrong use of uninitialized accessor.\n");
    return mesh_->element_accessor( elem_idx_ );
}

inline unsigned int Side::edge_idx() const {
    return element()->edge_idx(side_idx_);
}

inline Edge Side::edge() const {
    return mesh_->edge(edge_idx());
}

inline Boundary Side::cond() const {
    return mesh_->boundary(cond_idx());
}

inline unsigned int Side::cond_idx() const {
        if (element()->boundary_idx_ == nullptr) return Mesh::undef_idx;
        else return element()->boundary_idx_[side_idx_];
}



/*******************************************************************************
 * Boundary IMPLEMENTATION
 *******************************************************************************/
inline Boundary::Boundary()
: boundary_data_(nullptr)
{}

inline Boundary::Boundary(BoundaryData* boundary_data)
: boundary_data_(boundary_data)
{}

inline Edge Boundary::edge()
{
    ASSERT_DBG(is_valid());
    return boundary_data_->mesh_->edge(boundary_data_->edge_idx_);
}

inline ElementAccessor<3> Boundary::element_accessor()
{
    ASSERT_DBG(is_valid());
    return boundary_data_->mesh_->element_accessor(boundary_data_->bc_ele_idx_);
}

inline Region Boundary::region()
{ 
    return element_accessor().region();
}

inline Element * Boundary::element()
{
    ASSERT_DBG(is_valid());
    return &( boundary_data_->mesh_->element_vec_[boundary_data_->bc_ele_idx_] );
}