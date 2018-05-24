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

#include "mesh/side_impl.hh"
#include "output_mesh.hh"
#include "output_element.hh"
#include "mesh/mesh.h"
#include "mesh/ref_element.hh"
#include "mesh/long_idx.hh"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"
#include "la/distribution.hh"


namespace IT=Input::Type;

const IT::Record & OutputMeshBase::get_input_type() {
    return IT::Record("OutputMesh", "Parameters of the refined output mesh.")
        .declare_key("max_level", IT::Integer(1,20),IT::Default("3"),
            "Maximal level of refinement of the output mesh.")
        .declare_key("refine_by_error", IT::Bool(), IT::Default("false"),
            "Set true for using error_control_field. Set false for global uniform refinement to max_level.")
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
	mesh_type_(MeshType::orig),
    refine_by_error_(false),
    refinement_error_tolerance_(0.0)
{
}


OutputMeshBase::OutputMeshBase(Mesh &mesh, const Input::Record &in_rec)
: 
    input_record_(in_rec), 
    orig_mesh_(&mesh),
    max_level_(input_record_.val<int>("max_level")),
	mesh_type_(MeshType::orig),
    refine_by_error_(input_record_.val<bool>("refine_by_error")),
    refinement_error_tolerance_(input_record_.val<double>("refinement_error_tolerance"))
{
}

OutputMeshBase::~OutputMeshBase()
{
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
	elem_ids_ = std::make_shared< ElementDataCache<unsigned int> >("elements_ids", (unsigned int)1, 1, this->n_elements());
	node_ids_ = std::make_shared< ElementDataCache<unsigned int> >("node_ids", (unsigned int)1, 1, this->n_nodes());
	region_ids_ = std::make_shared< ElementDataCache<unsigned int> >("region_ids", (unsigned int)1, 1, this->n_elements());
	partitions_ = std::make_shared< ElementDataCache<int> >("partitions", (unsigned int)1, 1, this->n_elements());
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
			if (mesh_type_ == MeshType::orig) node_idx[0] = orig_mesh_->node_vector(node_list[j]).id();
			else node_idx[0] = node_list[j];
			node_ids_->store_value( node_list[j], node_idx );
		}
	}
}


bool OutputMeshBase::is_created()
{
	return (nodes_ && connectivity_ && offsets_);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////


OutputMesh::OutputMesh(Mesh  &mesh)
: OutputMeshBase(mesh)
{
}

OutputMesh::OutputMesh(Mesh &mesh, const Input::Record& in_rec)
: OutputMeshBase(mesh, in_rec)
{
}


OutputMesh::~OutputMesh()
{
}


void OutputMesh::create_mesh()
{
	ASSERT( !is_created() ).error("Multiple initialization of OutputMesh!\n");

	DebugOut() << "Create outputmesh identical to computational one.";

    const unsigned int n_elements = orig_mesh_->n_elements(),
                       n_nodes = orig_mesh_->n_nodes();

    nodes_ = std::make_shared<ElementDataCache<double>>("", (unsigned int)ElementDataCacheBase::N_VECTOR, 1, n_nodes);
    unsigned int coord_id = 0,  // coordinate id in vector
                 node_id = 0;   // node id
    auto &node_vec = *( nodes_->get_component_data(0).get() );
    FOR_NODES(orig_mesh_, node) {
        node->aux = node_id;   // store node index in the auxiliary variable

        // use store value
        node_vec[coord_id] = node->getX();  coord_id++;
        node_vec[coord_id] = node->getY();  coord_id++;
        node_vec[coord_id] = node->getZ();  coord_id++;
        node_id++;
    }
    
    orig_element_indices_ = std::make_shared<std::vector<unsigned int>>(n_elements);

    offsets_ = std::make_shared<ElementDataCache<unsigned int>>("offsets", (unsigned int)ElementDataCacheBase::N_SCALAR, 1, n_elements);
    Node* node;
    unsigned int ele_id = 0,
                 connect_id = 0,
                 offset = 0,    // offset of node indices of element in node vector
                 li;            // local node index
    auto &offset_vec = *( offsets_->get_component_data(0).get() );
    for (auto ele : orig_mesh_->bulk_elements_range()) {
        // increase offset by number of nodes of the simplicial element
        offset += ele->dim() + 1;
        offset_vec[ele_id] = offset;
        (*orig_element_indices_)[ele_id] = ele_id;
        ele_id++;
    }

    const unsigned int n_connectivities = offset_vec[offset_vec.size()-1];
    connectivity_ = std::make_shared<ElementDataCache<unsigned int>>("connectivity", (unsigned int)ElementDataCacheBase::N_SCALAR,
    		1, n_connectivities);
    auto &connect_vec = *( connectivity_->get_component_data(0).get() );
    for (auto ele : orig_mesh_->bulk_elements_range()) {
    	for (li=0; li<ele->n_nodes(); li++) {
            node = ele->node[li];
            connect_vec[connect_id] = node->aux;
            connect_id++;
        }
        
    }
}

void OutputMesh::create_refined_mesh()
{
    ASSERT(0).error("Not implemented yet.");
}


bool OutputMesh::refinement_criterion()
{
    ASSERT(0).error("Not implemented yet.");
    return false;
}


void OutputMesh::create_sub_mesh()
{
    ASSERT(0).error("Not implemented yet.");
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////

OutputMeshDiscontinuous::OutputMeshDiscontinuous(Mesh &mesh)
: OutputMeshBase(mesh)
{
}

OutputMeshDiscontinuous::OutputMeshDiscontinuous(Mesh &mesh, const Input::Record& in_rec)
: OutputMeshBase(mesh, in_rec)
{
}


OutputMeshDiscontinuous::~OutputMeshDiscontinuous()
{
}


void OutputMeshDiscontinuous::create_mesh()
{
	ASSERT( !is_created() ).error("Multiple initialization of OutputMesh!\n");

    ASSERT_DBG(orig_mesh_->n_nodes() > 0);   //continuous data already computed
    
    if (nodes_) return;          //already computed
    
    DebugOut() << "Create discontinuous outputmesh.";
    
    const unsigned int n_elements = orig_mesh_->n_elements();

    orig_element_indices_ = std::make_shared<std::vector<unsigned int>>(n_elements);
    offsets_ = std::make_shared<ElementDataCache<unsigned int>>("offsets", (unsigned int)ElementDataCacheBase::N_SCALAR, 1, n_elements);

    unsigned int ele_id = 0,
                 offset = 0,    // offset of node indices of element in node vector
                 coord_id = 0,  // coordinate id in vector
                 corner_id = 0, // corner index (discontinous node)
                 li = 0;        // local node index

    auto &offset_vec = *( offsets_->get_component_data(0).get() );
    for (auto ele : orig_mesh_->bulk_elements_range()) {
        // increase offset by number of nodes of the simplicial element
        offset += ele->dim() + 1;
        offset_vec[ele_id] = offset;
        (*orig_element_indices_)[ele_id] = ele_id;
        ele_id++;
    }

    // connectivity = for every element list the nodes => its length corresponds to discontinuous data
    const unsigned int n_corners = offset_vec[offset_vec.size()-1];

    nodes_ = std::make_shared<ElementDataCache<double>>("", (unsigned int)ElementDataCacheBase::N_VECTOR, 1, n_corners);
    connectivity_ = std::make_shared<ElementDataCache<unsigned int>>("connectivity", (unsigned int)ElementDataCacheBase::N_SCALAR,
    		1, n_corners);

    auto &node_vec = *( nodes_->get_component_data(0).get() );
    auto &conn_vec = *( connectivity_->get_component_data(0).get() );
    Node* node;
    for (auto ele : orig_mesh_->bulk_elements_range()) {
    	for (li=0; li<ele->n_nodes(); li++)
        {
            node = ele->node[li];
            node_vec[coord_id] = node->getX();  ++coord_id;
            node_vec[coord_id] = node->getY();  ++coord_id;
            node_vec[coord_id] = node->getZ();  ++coord_id;

            conn_vec[corner_id] = corner_id;
            corner_id++;
        }
    }

    mesh_type_ = MeshType::discont;
}


void OutputMeshDiscontinuous::create_refined_mesh()
{
    DebugOut() << "Create refined discontinuous outputmesh.\n";
    // initial guess of size: n_elements
    nodes_ = std::make_shared<ElementDataCache<double>>("",(unsigned int)ElementDataCacheBase::N_VECTOR,1,0);
    connectivity_ = std::make_shared<ElementDataCache<unsigned int>>("connectivity",(unsigned int)ElementDataCacheBase::N_SCALAR,1,0);
    offsets_ = std::make_shared<ElementDataCache<unsigned int>>("offsets",(unsigned int)ElementDataCacheBase::N_SCALAR,1,0);
    orig_element_indices_ = std::make_shared<std::vector<unsigned int>>();
    
    // index of last node added; set at the end of original ones
    unsigned int last_offset = 0;
    
    auto &node_vec = *( nodes_->get_component_data(0).get() );
    auto &conn_vec = *( connectivity_->get_component_data(0).get() );
    auto &offset_vec = *( offsets_->get_component_data(0).get() );
    
    node_vec.reserve(4*orig_mesh_->n_nodes());
    conn_vec.reserve(4*4*orig_mesh_->n_elements());
    offset_vec.reserve(4*orig_mesh_->n_elements());
    
//     DebugOut() << "start refinement\n";
    for (auto ele : orig_mesh_->bulk_elements_range()) {
        const unsigned int
            dim = ele->dim(),
            ele_idx = ele.index();
//         DebugOut() << "ele index " << ele_idx << "\n";
        
        AuxElement aux_ele;
        aux_ele.nodes.resize(ele->n_nodes());
        aux_ele.level = 0;
        
        Node* node; unsigned int li;
        for (li=0; li<ele->n_nodes(); li++) {
            node = ele->node[li];
            aux_ele.nodes[li] = node->point();
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
    
    mesh_type_ = MeshType::refined;
//     for(unsigned int i=0; i< nodes_->n_values; i++)
//     {
//         cout << i << "  ";
//         for(unsigned int k=0; k<spacedim; k++){
//             nodes_->print(cout, i*spacedim+k);
//             cout << " ";
//         }
//         cout << endl;
//     }
//     cout << "\n\n";
// //     nodes_->print_all(cout);
// //     cout << "\n\n";
//     connectivity_->print_all(cout);
//     cout << "\n\n";
//     offsets_->print_all(cout);
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

void OutputMeshDiscontinuous::create_sub_mesh()
{
	ASSERT( !is_created() ).error("Multiple initialization of OutputMesh!\n");

	DebugOut() << "Create output submesh containing only local elements.";

	ElementAccessor<3> ele;
	LongIdx *el_4_loc = orig_mesh_->get_el_4_loc();
	Distribution *el_ds = orig_mesh_->get_el_ds();
    const unsigned int n_local_elements = el_ds->lsize();

    orig_element_indices_ = std::make_shared<std::vector<unsigned int>>(n_local_elements);
    offsets_ = std::make_shared<ElementDataCache<unsigned int>>("offsets", (unsigned int)ElementDataCacheBase::N_SCALAR, 1, n_local_elements);

    unsigned int ele_id = 0,
                 offset = 0,    // offset of node indices of element in node vector
                 coord_id = 0,  // coordinate id in vector
                 corner_id = 0, // corner index (discontinous node)
                 li = 0;        // local node index
    auto &offset_vec = *( offsets_->get_component_data(0).get() );
	for (unsigned int loc_el = 0; loc_el < n_local_elements; loc_el++) {
		ele = orig_mesh_->element_accessor( el_4_loc[loc_el] );
        // increase offset by number of nodes of the simplicial element
        offset += ele->dim() + 1;
        offset_vec[ele_id] = offset;
        (*orig_element_indices_)[ele_id] = el_4_loc[loc_el];
        ele_id++;
	}

    // connectivity = for every element list the nodes => its length corresponds to discontinuous data
    const unsigned int n_corners = offset_vec[offset_vec.size()-1];

    nodes_ = std::make_shared<ElementDataCache<double>>("", (unsigned int)ElementDataCacheBase::N_VECTOR, 1, n_corners);
    connectivity_ = std::make_shared<ElementDataCache<unsigned int>>("connectivity", (unsigned int)ElementDataCacheBase::N_SCALAR,
    		1, n_corners);

    auto &node_vec = *( nodes_->get_component_data(0).get() );
    auto &conn_vec = *( connectivity_->get_component_data(0).get() );
	Node* node;
	for (unsigned int loc_el = 0; loc_el < n_local_elements; loc_el++) {
		ele = orig_mesh_->element_accessor( el_4_loc[loc_el] );
		for (li=0; li<ele->n_nodes(); li++)
        {
            node = ele->node[li];
            node_vec[coord_id] = node->getX();  ++coord_id;
            node_vec[coord_id] = node->getY();  ++coord_id;
            node_vec[coord_id] = node->getZ();  ++coord_id;

            conn_vec[corner_id] = corner_id;
            corner_id++;
        }
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
    std::vector< Space<spacedim>::Point > point_list;
    point_list.push_back(centre);
    point_list.insert(point_list.end(), ele.nodes.begin(), ele.nodes.end());
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
