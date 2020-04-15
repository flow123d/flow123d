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
 * @file    output_mesh.cc
 * @brief   Classes for auxiliary output mesh.
 */

#include "system/index_types.hh"
#include "output_mesh.hh"
#include "output_element.hh"
#include "mesh/mesh.h"
#include "mesh/ref_element.hh"
#include "mesh/accessors.hh"
#include "mesh/node_accessor.hh"
#include "mesh/range_wrapper.hh"
#include "la/distribution.hh"


namespace IT=Input::Type;

const IT::Record & OutputMeshBase::get_input_type() {
    return IT::Record("OutputMesh", "Parameters of the refined output mesh. [Not impemented]")
        .declare_key("max_level", IT::Integer(1,20),IT::Default("3"),
            "Maximal level of refinement of the output mesh.")
        .declare_key("refine_by_error", IT::Bool(), IT::Default("false"),
            "Set true for using ``error_control_field``. Set false for global uniform refinement to max_level.")
        .declare_key("error_control_field",IT::String(), IT::Default::optional(),
            "Name of an output field, according to which the output mesh will be refined. The field must be a SCALAR one.")
        .declare_key("refinement_error_tolerance",IT::Double(0.0), IT::Default("0.01"),
            "Tolerance for element refinement by error. If tolerance is reached, refinement is stopped."
            "Relative difference between error control field and its linear approximation on element is computed"
            "and compared with tolerance.")
        .close();
}

OutputMeshBase::OutputMeshBase(Mesh &mesh)
: 
	orig_mesh_(&mesh),
    max_level_(0),
    refine_by_error_(false),
    refinement_error_tolerance_(0.0),
	el_ds_(nullptr),
	node_ds_(nullptr)
{
}


OutputMeshBase::OutputMeshBase(Mesh &mesh, const Input::Record &in_rec)
: 
    input_record_(in_rec), 
    orig_mesh_(&mesh),
    max_level_(input_record_.val<int>("max_level")),
    refine_by_error_(input_record_.val<bool>("refine_by_error")),
    refinement_error_tolerance_(input_record_.val<double>("refinement_error_tolerance")),
	el_ds_(nullptr),
	node_ds_(nullptr)
{
}

OutputMeshBase::~OutputMeshBase()
{
	// Refined mesh creates own special distributions and local to global maps and needs destroy these objects.
    if ( (el_ds_!=nullptr) && (mesh_type_ == MeshType::refined)) {
    	delete[] el_4_loc_;
    	delete[] node_4_loc_;
    	delete el_ds_;
    	delete node_ds_;
    }
}

OutputElementIterator OutputMeshBase::begin()
{
    ASSERT_PTR_DBG(offsets_);
//     ASSERT_DBG(offsets_->n_values() > 0);
    return OutputElementIterator(OutputElement(0, shared_from_this()));
}

OutputElementIterator OutputMeshBase::end()
{
    ASSERT_PTR_DBG(offsets_);
//     ASSERT_DBG(offsets_->n_values() > 0);
    return OutputElementIterator(OutputElement(offsets_->n_values(), shared_from_this()));
}


void OutputMeshBase::set_error_control_field(ErrorControlFieldFunc error_control_field_func)
{
    error_control_field_func_ = error_control_field_func;
}

unsigned int OutputMeshBase::n_elements()
{
    ASSERT_PTR(offsets_);
    return offsets_->n_values();
}

unsigned int OutputMeshBase::n_nodes()
{
    ASSERT_PTR(nodes_);
    return nodes_->n_values();
}

void OutputMeshBase::create_id_caches()
{
	unsigned int elm_idx[1];
	unsigned int node_idx[1];
	unsigned int region_idx[1];
	int partition[1];
	elem_ids_ = std::make_shared< ElementDataCache<unsigned int> >("elements_ids", (unsigned int)1, this->n_elements());
	node_ids_ = std::make_shared< ElementDataCache<unsigned int> >("node_ids", (unsigned int)1, this->n_nodes());
	region_ids_ = std::make_shared< ElementDataCache<unsigned int> >("region_ids", (unsigned int)1, this->n_elements());
	partitions_ = std::make_shared< ElementDataCache<int> >("partitions", (unsigned int)1, this->n_elements());
	OutputElementIterator it = this->begin();
	for (unsigned int i = 0; i < this->n_elements(); ++i, ++it) {
		if (mesh_type_ == MeshType::orig) elm_idx[0] = orig_mesh_->find_elem_id(it->idx());
		else elm_idx[0] = it->idx();
		elem_ids_->store_value( i, elm_idx );

		region_idx[0] = orig_mesh_->element_accessor(it->idx()).region().id();
		region_ids_->store_value( i, region_idx );

		partition[0] = orig_mesh_->element_accessor(it->idx())->pid();
		partitions_->store_value( i, partition );

		std::vector< unsigned int > node_list = it->node_list();
		for (unsigned int j = 0; j < it->n_nodes(); ++j) {
			if (mesh_type_ == MeshType::orig) node_idx[0] = orig_mesh_->find_node_id(node_list[j]);
			else node_idx[0] = node_list[j];
			node_ids_->store_value( node_list[j], node_idx );
		}
	}
}


bool OutputMeshBase::is_created()
{
	return (nodes_ && connectivity_ && offsets_);
}


void OutputMeshBase::create_sub_mesh()
{
	ASSERT( !is_created() ).error("Multiple initialization of OutputMesh!\n");

	DebugOut() << "Create output submesh containing only local elements.";

    unsigned int ele_id = 0,
                 offset = 0,    // offset of node indices of element in node vector
                 coord_id = 0,  // coordinate id in node vector
                 conn_id = 0;   // index to connectivity vector
    ElementAccessor<3> elm;

    el_4_loc_ = orig_mesh_->get_el_4_loc();
    el_ds_ = orig_mesh_->get_el_ds();
    node_4_loc_ = orig_mesh_->get_node_4_loc();
    node_ds_ = orig_mesh_->get_node_ds();
    n_local_nodes_ = orig_mesh_->n_local_nodes();

    const unsigned int n_local_elements = el_ds_->lsize();
    unsigned int n_nodes = node_ds_->end( node_ds_->np()-1 );
    std::vector<unsigned int> local_nodes_map(n_nodes, Mesh::undef_idx); // map global to local ids of nodes
    for (unsigned int i=0; i<n_local_nodes_; ++i) local_nodes_map[ node_4_loc_[i] ] = i;

    orig_element_indices_ = std::make_shared<std::vector<unsigned int>>(n_local_elements);
    offsets_ = std::make_shared<ElementDataCache<unsigned int>>("offsets", (unsigned int)ElementDataCacheBase::N_SCALAR, n_local_elements);
    auto &offset_vec = *( offsets_->get_component_data(0).get() );

    for (unsigned int loc_el = 0; loc_el < n_local_elements; loc_el++) {
        elm = orig_mesh_->element_accessor( el_4_loc_[loc_el] );
        // increase offset by number of nodes of the simplicial element
        offset += elm->dim() + 1;
        offset_vec[ele_id] = offset;
        (*orig_element_indices_)[ele_id] = el_4_loc_[loc_el];
        ele_id++;
    }

    connectivity_ = std::make_shared<ElementDataCache<unsigned int>>("connectivity", (unsigned int)ElementDataCacheBase::N_SCALAR,
            offset_vec[offset_vec.size()-1]);
    auto &connectivity_vec = *( connectivity_->get_component_data(0).get() );
    for (unsigned int loc_el = 0; loc_el < n_local_elements; loc_el++) {
        elm = orig_mesh_->element_accessor( el_4_loc_[loc_el] );
        for (unsigned int li=0; li<elm->n_nodes(); li++) {
        	ASSERT_DBG(local_nodes_map[ elm.node(li).idx() ] != Mesh::undef_idx)(elm.node(li).idx()).error("Undefined global to local node index!");
        	connectivity_vec[conn_id++] = local_nodes_map[ elm.node(li).idx() ];
        }
    }

    // set coords of nodes
    nodes_ = std::make_shared<ElementDataCache<double>>("", (unsigned int)ElementDataCacheBase::N_VECTOR, n_local_nodes_);
    auto &node_vec = *( nodes_->get_component_data(0) );
    for(unsigned int i_node=0; i_node<local_nodes_map.size(); ++i_node) {
        if (local_nodes_map[i_node]==Mesh::undef_idx) continue; // skip element if it is not local
        auto node = *orig_mesh_->node(i_node);
        coord_id = 3*local_nodes_map[i_node]; // id of first coordinates in node_vec
        node_vec[coord_id++] = node[0];
        node_vec[coord_id++] = node[1];
        node_vec[coord_id] = node[2];
    }
}



void OutputMeshBase::make_serial_master_mesh()
{
    std::shared_ptr<ElementDataCache<unsigned int>> global_offsets; // needs for creating serial nodes and connectivity caches on zero process
    auto elems_n_nodes = get_elems_n_nodes(); // collects number of nodes on each elements (for fill master_mesh_->offsets_)
    int rank = el_ds_->myp();

    if (rank==0) {
    	// create serial output mesh, fill offsets cache and orig_element_indices vector
    	unsigned int n_elems = el_ds_->end( el_ds_->np()-1 );
    	master_mesh_ = this->construct_mesh();
    	master_mesh_->orig_element_indices_ = std::make_shared<std::vector<unsigned int>>(n_elems);
    	master_mesh_->offsets_ = std::make_shared<ElementDataCache<unsigned int>>("offsets", ElementDataCacheBase::N_SCALAR, n_elems);
        auto &offsets_vec = *( master_mesh_->offsets_->get_component_data(0).get() );
        auto &elems_n_nodes_vec = *( elems_n_nodes->get_component_data(0).get() );
        unsigned int offset=0;
        for (unsigned int i=0; i<n_elems; ++i) {
            offset += elems_n_nodes_vec[i];
            offsets_vec[i] = offset;
            (*master_mesh_->orig_element_indices_)[i] = i;
        }
        global_offsets = master_mesh_->offsets_;
    }

    // collects serial caches
    std::shared_ptr<ElementDataCache<double>> serial_nodes_cache = make_serial_nodes_cache(global_offsets);
    std::shared_ptr<ElementDataCache<unsigned int>> serial_connectivity_cache = make_serial_connectivity_cache(global_offsets);

    if (rank==0) {
        // set serial output mesh caches
    	master_mesh_->connectivity_ = serial_connectivity_cache;
    	master_mesh_->nodes_ = serial_nodes_cache;

    	master_mesh_->mesh_type_ = this->mesh_type_;
    }
}


std::shared_ptr<ElementDataCache<unsigned int>> OutputMeshBase::get_elems_n_nodes()
{
	// Compute (locally) number of nodes of each elements
	ElementDataCache<unsigned int> local_elems_n_nodes("elems_n_nodes", ElementDataCacheBase::N_SCALAR, offsets_->n_values());
	auto &local_elems_n_nodes_vec = *( local_elems_n_nodes.get_component_data(0).get() );
	auto &offset_vec = *( offsets_->get_component_data(0).get() );
	for (unsigned int i=offset_vec.size()-1; i>0; --i) local_elems_n_nodes_vec[i] = offset_vec[i] - offset_vec[i-1];
	local_elems_n_nodes_vec[0] = offset_vec[0];

	// Collect data, set on zero process
	std::shared_ptr<ElementDataCache<unsigned int>> global_elems_n_nodes;
	auto gather_cache = local_elems_n_nodes.gather(el_ds_, el_4_loc_);
	if (el_ds_->myp()==0) global_elems_n_nodes = std::dynamic_pointer_cast< ElementDataCache<unsigned int> >(gather_cache);
	return global_elems_n_nodes;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////


OutputMesh::OutputMesh(Mesh  &mesh)
: OutputMeshBase(mesh)
{
    this->mesh_type_ = MeshType::orig;
}

OutputMesh::OutputMesh(Mesh &mesh, const Input::Record& in_rec)
: OutputMeshBase(mesh, in_rec)
{
    this->mesh_type_ = MeshType::orig;
}


OutputMesh::~OutputMesh()
{
}


void OutputMesh::create_refined_sub_mesh()
{
    ASSERT(0).error("Not implemented yet.");
}

bool OutputMesh::refinement_criterion()
{
    ASSERT(0).error("Not implemented yet.");
    return false;
}


std::shared_ptr<OutputMeshBase> OutputMesh::construct_mesh()
{
    return std::make_shared<OutputMesh>(*orig_mesh_);
}


std::shared_ptr<ElementDataCache<double>> OutputMesh::make_serial_nodes_cache(FMT_UNUSED std::shared_ptr<ElementDataCache<unsigned int>> global_offsets)
{
	std::shared_ptr<ElementDataCache<double>> serial_nodes_cache;

    // collects nodes_ data (coordinates)
    auto serial_nodes = nodes_->gather(node_ds_, node_4_loc_);

    if (el_ds_->myp()==0) serial_nodes_cache = std::dynamic_pointer_cast< ElementDataCache<double> >(serial_nodes);
    return serial_nodes_cache;
}


std::shared_ptr<ElementDataCache<unsigned int>> OutputMesh::make_serial_connectivity_cache(std::shared_ptr<ElementDataCache<unsigned int>> global_offsets)
{
	std::shared_ptr<ElementDataCache<unsigned int>> serial_connectivity_cache;

    // re-number connectivity indices from local to global
    auto &conn_vec = *( connectivity_->get_component_data(0).get() );
    ElementDataCache<unsigned int> global_conn("connectivity", (unsigned int)1, conn_vec.size()); // holds global indices of nodes
    auto &global_conn_vec = *( global_conn.get_component_data(0).get() );
    for(unsigned int i=0; i<conn_vec.size(); i++) {
        global_conn_vec[i] = node_4_loc_[ conn_vec[i] ];
    }

    // collects global connectivities
    auto &local_offset_vec = *( offsets_->get_component_data(0).get() );
    auto global_fix_size_conn = global_conn.element_node_cache_fixed_size(local_offset_vec);
    auto collective_conn = global_fix_size_conn->gather(el_ds_, el_4_loc_);

    if (el_ds_->myp()==0) {
    	auto &offset_vec = *( global_offsets->get_component_data(0).get() );
    	serial_connectivity_cache = std::dynamic_pointer_cast< ElementDataCache<unsigned int> >( collective_conn->element_node_cache_optimize_size(offset_vec) );
    }
    return serial_connectivity_cache;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

OutputMeshDiscontinuous::OutputMeshDiscontinuous(Mesh &mesh)
: OutputMeshBase(mesh)
{
    this->mesh_type_ = MeshType::discont;
}

OutputMeshDiscontinuous::OutputMeshDiscontinuous(Mesh &mesh, const Input::Record& in_rec)
: OutputMeshBase(mesh, in_rec)
{
    this->mesh_type_ = MeshType::discont;
}


OutputMeshDiscontinuous::~OutputMeshDiscontinuous()
{
}


template<int dim>
void OutputMeshDiscontinuous::refine_aux_element(const OutputMeshDiscontinuous::AuxElement& aux_element,
                                                 std::vector< OutputMeshDiscontinuous::AuxElement >& refinement,
                                                 const ElementAccessor<spacedim> &ele_acc)
{
    static const unsigned int n_subelements = 1 << dim;  //2^dim
    
// The refinement of elements for the output mesh is done using edge splitting
// technique (so called red refinement). Since we use this only for better output
// visualization of non-polynomial solutions, we do not care for existence of hanging
// nodes.
// In 2D case, it is straightforward process: find the midpoints of all sides,
// connect them and generate 4 triangles. These triangles are congruent and have
// equal surface areas.
// On the other hand, the 3D case is more complicated. After splitting the
// edges, we obtain 4 tetrahedra at the vertices of the original one. The octahedron
// that remains in the middle can be subdivided according to one of its three
// diagonals. Only the choice of the shortest octahedron diagonal leads to a regular
// tetrahedra decomposition. This algorithm originally comes from Bey.
//  Bey's algorithm (red refinement of tetrahedron):
// p.29 https://www5.in.tum.de/pub/Joshi2016_Thesis.pdf
// p.108 http://www.bcamath.org/documentos_public/archivos/publicaciones/sergey_book.pdf
// https://www.math.uci.edu/~chenlong/iFEM/doc/html/uniformrefine3doc.html#1
// J. Bey. Simplicial grid refinement: on Freudenthal's algorithm and the optimal number of congruence classes.
//    Numer. Math. 85(1):1--29, 2000. p11 Algorithm: RedRefinement3D.
// p.4 http://www.vis.uni-stuttgart.de/uploads/tx_vispublications/vis97-grosso.pdf

    // connectivity of refined element
    // these arrays are hardwired to the current reference element
    static const std::vector<unsigned int> conn[] = {
        {}, //0D
        
        //1D:
        // 0,1 original nodes, 2 is in the middle
        // get 2 elements
        {0, 2,
         2, 1},
        
        //2D:
        // 0,1,2 original nodes
        // 3,4,5 nodes are in the middle of sides 0,1,2 in the respective order
        // get 4 elements
        {0, 3, 4,
         3, 1, 5,
         4, 5, 2,
         3, 5, 4},
        
        //3D:
        // 0,1,2,3 original nodes
        // 4,5,6,7,8,9 are nodes in the middle of edges 0,1,2,3,4,5 in the respective order 
        // first 4 tetrahedra are at the original nodes
        // next 4 tetrahedra are from the inner octahedron - 4 choices according to the diagonal
        {1, 7, 4, 8,
         7, 2, 5, 9,
         4, 5, 0, 6,
         8, 9, 6, 3,
         7, 4, 8, 9, // 4-9 octahedron diagonal
         7, 4, 5, 9,
         4, 8, 9, 6,
         4, 5, 9, 6},
        
        {1, 7, 4, 8,
         7, 2, 5, 9,
         4, 5, 0, 6,
         8, 9, 6, 3,
         8, 5, 4, 6, // 5-8 octahedron diagonal
         5, 8, 9, 6,
         5, 8, 4, 7,
         8, 5, 9, 7},
         
        {1, 7, 4, 8,
         7, 2, 5, 9,
         4, 5, 0, 6,
         8, 9, 6, 3,
         6, 8, 7, 9, // 6-7 octahedron diagonal
         8, 6, 7, 4,
         6, 5, 7, 4,
         5, 6, 7, 9}
    };
//     DBGMSG("level = %d, %d\n", aux_element.level, max_refinement_level_);
 
    ASSERT_DBG(dim == aux_element.nodes.size()-1);
    
    // if not refining any further, push into final vector
    if( ! refinement_criterion(aux_element, ele_acc) ) {
        refinement.push_back(aux_element);
        return;
    }
    
    std::vector<AuxElement> subelements(n_subelements);
    
    const unsigned int n_old_nodes = RefElement<dim>::n_nodes,
                       n_new_nodes = RefElement<dim>::n_lines; // new points are in the center of lines
    
    // auxiliary vectors
    std::vector<Space<spacedim>::Point> nodes = aux_element.nodes;
    
    nodes.reserve(n_old_nodes+n_new_nodes);

    // create new points in the element
    for(unsigned int e=0; e < n_new_nodes; e++)
    {
        Space<spacedim>::Point p = nodes[RefElement<dim>::interact(Interaction<0,1>(e))[0]]
                                  +nodes[RefElement<dim>::interact(Interaction<0,1>(e))[1]];
        nodes.push_back( p / 2.0);
        //nodes.back().print();
    }
    
    unsigned int diagonal = 0;
    // find shortest diagonal: [0]:4-9, [1]:5-8 or [2]:6-7
    if(dim == 3){
        double min_diagonal = arma::norm(nodes[4]-nodes[9],2);
        double d = arma::norm(nodes[5]-nodes[8],2);
        if(d < min_diagonal){
            min_diagonal = d;
            diagonal = 1;
        }
        d = arma::norm(nodes[6]-nodes[7],2);
        if(d < min_diagonal){
            min_diagonal = d;
            diagonal = 2;
        }
    }
    
    for(unsigned int i=0; i < n_subelements; i++)
    {
        AuxElement& sub_ele = subelements[i];
        sub_ele.nodes.resize(n_old_nodes);
        sub_ele.level = aux_element.level+1;
        
        // over nodes
        for(unsigned int j=0; j < n_old_nodes; j++)
        {
            unsigned int conn_id = (n_old_nodes)*i + j;
            sub_ele.nodes[j] = nodes[conn[dim+diagonal][conn_id]];
        }
        refine_aux_element<dim>(sub_ele, refinement, ele_acc);
    }
}



template void OutputMeshDiscontinuous::refine_aux_element<1>(const OutputMeshDiscontinuous::AuxElement&,std::vector< OutputMeshDiscontinuous::AuxElement >&, const ElementAccessor<spacedim> &);
template void OutputMeshDiscontinuous::refine_aux_element<2>(const OutputMeshDiscontinuous::AuxElement&,std::vector< OutputMeshDiscontinuous::AuxElement >&, const ElementAccessor<spacedim> &);
template void OutputMeshDiscontinuous::refine_aux_element<3>(const OutputMeshDiscontinuous::AuxElement&,std::vector< OutputMeshDiscontinuous::AuxElement >&, const ElementAccessor<spacedim> &);


bool OutputMeshDiscontinuous::refinement_criterion(const AuxElement& aux_ele,
                                                   const ElementAccessor<spacedim> &ele_acc)
{
    // check refinement criteria:
    
    //first check max. level
    bool refine = refinement_criterion_uniform(aux_ele);
    
    //if max. level not reached and refinement by error is set
    if(refine && refine_by_error_)
    {
        // compute centre of aux element
        Space<spacedim>::Point centre({0,0,0});
        for(auto& v : aux_ele.nodes ) centre += v;
        centre = centre/aux_ele.nodes.size();
        return refinement_criterion_error(aux_ele, centre, ele_acc);
    }
    
    return refine;
}

bool OutputMeshDiscontinuous::refinement_criterion_uniform(const OutputMeshDiscontinuous::AuxElement& ele)
{
    return (ele.level < max_level_);
}

bool OutputMeshDiscontinuous::refinement_criterion_error(const OutputMeshDiscontinuous::AuxElement& ele,
                                            const Space<spacedim>::Point &centre,
                                            const ElementAccessor<spacedim> &ele_acc
                                           )
{
    ASSERT_DBG(error_control_field_func_).error("Error control field not set!");

    // evaluate at nodes and center in a single call
    std::vector<double> val_list(ele.nodes.size()+1);
    Armor::array point_list(spacedim,1,1+ele.nodes.size());
    point_list.set(0) = centre;
    unsigned int i=0;
    for (auto node : ele.nodes) point_list.set(++i) = node;
    error_control_field_func_(point_list, ele_acc, val_list);

    //TODO: compute L1 or L2 error using standard quadrature
    
    //compare average value at nodes with value at center
    
    double average_val = 0.0;
    for(unsigned int i=1; i<ele.nodes.size()+1; ++i)//(double& v: nodes_val)
        average_val += val_list[i];
    average_val = average_val / ele.nodes.size();
    
    double diff = std::abs((average_val - val_list[0])/val_list[0]);
//     DebugOut().fmt("diff: {}  {}  {}\n", diff, average_val, val_list[0]);
    return ( diff > refinement_error_tolerance_);

}


std::shared_ptr<OutputMeshBase> OutputMeshDiscontinuous::construct_mesh()
{
    return std::make_shared<OutputMeshDiscontinuous>(*orig_mesh_);
}


std::shared_ptr<ElementDataCache<double>> OutputMeshDiscontinuous::make_serial_nodes_cache(std::shared_ptr<ElementDataCache<unsigned int>> global_offsets)
{
	std::shared_ptr<ElementDataCache<double>> serial_nodes_cache;

    // Create helper cache of discontinuous node data ordering by elements
    std::shared_ptr< ElementDataCache<double> > discont_node_cache = std::make_shared<ElementDataCache<double>>("",
                ElementDataCacheBase::N_VECTOR, this->connectivity_->n_values());
    auto &discont_node_vec = *( discont_node_cache->get_component_data(0).get() );
    auto &local_nodes_vec = *( this->nodes_->get_component_data(0).get() );
    auto &local_conn_vec = *( this->connectivity_->get_component_data(0).get() );
    auto &local_offset_vec = *( this->offsets_->get_component_data(0).get() );
    unsigned int i_old, i_new;
    for (unsigned int i_conn=0; i_conn<this->connectivity_->n_values(); ++i_conn) {
    	i_old = local_conn_vec[i_conn] * ElementDataCacheBase::N_VECTOR;
    	i_new = i_conn * ElementDataCacheBase::N_VECTOR;
        for(unsigned int i = 0; i < ElementDataCacheBase::N_VECTOR; i++) {
        	discont_node_vec[i_new+i] = local_nodes_vec[i_old+i];
        }
    }
    // Collects node data
    auto fix_size_node_cache = discont_node_cache->element_node_cache_fixed_size(local_offset_vec);
    auto collect_fix_size_node_cache = fix_size_node_cache->gather(el_ds_, el_4_loc_);

    if (el_ds_->myp()==0) {
    	auto &offset_vec = *( global_offsets->get_component_data(0).get() );
        serial_nodes_cache = std::dynamic_pointer_cast< ElementDataCache<double> >(collect_fix_size_node_cache->element_node_cache_optimize_size(offset_vec));
    }
    return serial_nodes_cache;
}


std::shared_ptr<ElementDataCache<unsigned int>> OutputMeshDiscontinuous::make_serial_connectivity_cache(std::shared_ptr<ElementDataCache<unsigned int>> global_offsets)
{
	std::shared_ptr<ElementDataCache<unsigned int>> serial_connectivity_cache;

    if (el_ds_->myp()==0) {
    	auto &offset_vec = *( global_offsets->get_component_data(0).get() );
    	serial_connectivity_cache = std::make_shared<ElementDataCache<unsigned int>>("connectivity", (unsigned int)ElementDataCacheBase::N_SCALAR,
                offset_vec[offset_vec.size()-1]);
        auto &conn_vec = *( serial_connectivity_cache->get_component_data(0).get() );
        for (unsigned int i=0; i<conn_vec.size(); ++i) conn_vec[i] = i;
    }
    return serial_connectivity_cache;
}


void OutputMeshDiscontinuous::create_refined_sub_mesh()
{
    ASSERT( !is_created() ).error("Multiple initialization of OutputMesh!\n");

    DebugOut() << "Create refined discontinuous submesh containing only local elements.";
    // initial guess of size: n_elements
    nodes_ = std::make_shared<ElementDataCache<double>>("",(unsigned int)ElementDataCacheBase::N_VECTOR,0);
    connectivity_ = std::make_shared<ElementDataCache<unsigned int>>("connectivity",(unsigned int)ElementDataCacheBase::N_SCALAR,0);
    offsets_ = std::make_shared<ElementDataCache<unsigned int>>("offsets",(unsigned int)ElementDataCacheBase::N_SCALAR,0);
    orig_element_indices_ = std::make_shared<std::vector<unsigned int>>();

    // index of last node added; set at the end of original ones
    unsigned int last_offset = 0;

    auto &node_vec = *( nodes_->get_component_data(0).get() );
    auto &conn_vec = *( connectivity_->get_component_data(0).get() );
    auto &offset_vec = *( offsets_->get_component_data(0).get() );

    node_vec.reserve(4*orig_mesh_->n_nodes());
    conn_vec.reserve(4*4*orig_mesh_->n_elements());
    offset_vec.reserve(4*orig_mesh_->n_elements());

    LongIdx *el_4_loc = orig_mesh_->get_el_4_loc();
    const unsigned int n_local_elements = orig_mesh_->get_el_ds()->lsize();

    for (unsigned int loc_el = 0; loc_el < n_local_elements; loc_el++) {
    	auto ele = orig_mesh_->element_accessor( el_4_loc[loc_el] );
    	const unsigned int
            dim = ele->dim(),
            ele_idx = ele.idx();

        AuxElement aux_ele;
        aux_ele.nodes.resize(ele->n_nodes());
        aux_ele.level = 0;

        unsigned int li;
        for (li=0; li<ele->n_nodes(); li++) {
            aux_ele.nodes[li] = *ele.node(li);
        }

        std::vector<AuxElement> refinement;

        switch(dim){
            case 1: this->refine_aux_element<1>(aux_ele, refinement, ele); break;
            case 2: this->refine_aux_element<2>(aux_ele, refinement, ele); break;
            case 3: this->refine_aux_element<3>(aux_ele, refinement, ele); break;
            default: ASSERT(0 < dim && dim < 4);
        }

        //skip unrefined element
//         if(refinement.size() < 2) continue;
        unsigned int node_offset = node_vec.size(),
                     con_offset = conn_vec.size();
        node_vec.resize(node_vec.size() + (refinement.size() * (dim+1))*spacedim);
        conn_vec.resize(conn_vec.size() + refinement.size()*(dim+1));
//         orig_element_indices_->resize(orig_element_indices_->size() + refinement.size()*(dim+1));

//         DebugOut() << "ref size = " << refinement.size() << "\n";
        //gather coords and connectivity (in a continous way inside element)
        for(unsigned int i=0; i < refinement.size(); i++)
        {
            last_offset += dim+1;
            offset_vec.push_back(last_offset);
            (*orig_element_indices_).push_back(ele_idx);
            for(unsigned int j=0; j < dim+1; j++)
            {
                unsigned int con = i*(dim+1) + j;
                conn_vec[con_offset + con] = con_offset + con;

                for(unsigned int k=0; k < spacedim; k++) {
                    node_vec[node_offset + con*spacedim + k] = refinement[i].nodes[j][k];
                }
            }
        }
    }

    conn_vec.shrink_to_fit();
    node_vec.shrink_to_fit();
    offset_vec.shrink_to_fit();

    connectivity_->set_n_values(conn_vec.size());
    nodes_->set_n_values(node_vec.size() / spacedim);
    offsets_->set_n_values(offset_vec.size());

    // Create special distributions and arrays of local to global indexes of refined mesh
	el_ds_ = new Distribution(offset_vec.size(), PETSC_COMM_WORLD);
	node_ds_ = new Distribution(offset_vec[offset_vec.size()-1], PETSC_COMM_WORLD);
	n_local_nodes_ = node_ds_->lsize();
	el_4_loc_ = new LongIdx [ el_ds_->lsize() ];
	LongIdx global_el_idx = el_ds_->begin();
	for (unsigned int i=0; i<el_ds_->lsize(); ++i, ++global_el_idx) {
	    el_4_loc_[i] = global_el_idx;
	}
	node_4_loc_ = new LongIdx [ node_ds_->lsize() ];
	LongIdx global_node_idx = node_ds_->begin();
	for (unsigned int i=0; i<node_ds_->lsize(); ++i, ++global_node_idx) {
		node_4_loc_[i] = global_node_idx;
	}

	mesh_type_ = MeshType::refined;
}


void OutputMeshDiscontinuous::make_parallel_master_mesh()
{
    master_mesh_ = this->construct_mesh();
	master_mesh_->offsets_ = this->offsets_;
    master_mesh_->orig_element_indices_ = this->orig_element_indices_;

    auto &conn_vec = *( this->connectivity_->get_component_data(0).get() );
    master_mesh_->connectivity_ = std::make_shared<ElementDataCache<unsigned int>>("connectivity", ElementDataCacheBase::N_SCALAR,
            conn_vec.size());
    auto &master_conn_vec = *( master_mesh_->connectivity_->get_component_data(0).get() );
    for (unsigned int i=0; i<master_conn_vec.size(); ++i) master_conn_vec[i] = i;

    master_mesh_->nodes_ = std::make_shared<ElementDataCache<double>>("", ElementDataCacheBase::N_VECTOR, conn_vec.size());
    auto &node_vec = *( this->nodes_->get_component_data(0).get() );
    auto &master_node_vec = *( master_mesh_->nodes_->get_component_data(0).get() );
    unsigned int i_own, i_master, j;
    for (unsigned int i=0; i<conn_vec.size(); ++i) {
    	i_own = conn_vec[i]*ElementDataCacheBase::N_VECTOR;
    	i_master = i*ElementDataCacheBase::N_VECTOR;
    	for (j=0; j<ElementDataCacheBase::N_VECTOR; ++j) {
    		master_node_vec[i_master+j] = node_vec[i_own+j];
    	}
    }

	master_mesh_->mesh_type_ = this->mesh_type_;
}
