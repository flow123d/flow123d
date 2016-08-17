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

#include "output_mesh.hh"
#include "output_element.hh"
#include "mesh/mesh.h"
#include "mesh/ref_element.hh"
#include "fields/field.hh"
#include "fields/field_set.hh"

namespace IT=Input::Type;

const IT::Record & OutputMeshBase::get_input_type() {
    return IT::Record("OutputStream", "Parameters of output.")
        .declare_key("max_level", IT::Integer(1,20),IT::Default("3"),
            "Maximal level of refinement of the output mesh.")
        .declare_key("refine_by_error", IT::Bool(), IT::Default("false"),
            "Set true for using error_control_field. Set false for global uniform refinement to max_level.")
        .declare_key("error_control_field",IT::String(), IT::Default::optional(),
            "Name of an output field, according to which the output mesh will be refined. The field must be a SCALAR one.")
        .declare_key("refinement_error_tolerance",IT::Double(0.0), IT::Default("0.01"),
            "Tolerance for refinement by error.")
        .close();
}

OutputMeshBase::OutputMeshBase(Mesh &mesh)
: 
    nodes_ (std::make_shared<MeshData<double>>("", OutputDataBase::N_VECTOR)),
    connectivity_(std::make_shared<MeshData<unsigned int>>("connectivity")),
    offsets_(std::make_shared<MeshData<unsigned int>>("offsets")),
    orig_mesh_(&mesh),
    max_level_(0),
    is_refined_(false),
    refine_by_error_(false),
    refinement_error_tolerance_(0.0)
{
}


OutputMeshBase::OutputMeshBase(Mesh &mesh, const Input::Record &in_rec)
: 
    nodes_ (std::make_shared<MeshData<double>>("", OutputDataBase::N_VECTOR)),
    connectivity_(std::make_shared<MeshData<unsigned int>>("connectivity")),
    offsets_(std::make_shared<MeshData<unsigned int>>("offsets")),
    input_record_(in_rec), 
    orig_mesh_(&mesh),
    max_level_(input_record_.val<int>("max_level")),
    is_refined_(false),
    refine_by_error_(input_record_.val<bool>("refine_by_error")),
    refinement_error_tolerance_(input_record_.val<double>("refinement_error_tolerance"))
{
}

OutputMeshBase::~OutputMeshBase()
{
}

OutputElementIterator OutputMeshBase::begin()
{
    ASSERT_PTR(offsets_);
    ASSERT_DBG(offsets_->n_values > 0);
    return OutputElementIterator(OutputElement(0, shared_from_this()));
}

OutputElementIterator OutputMeshBase::end()
{
    ASSERT_PTR_DBG(offsets_);
    ASSERT_DBG(offsets_->n_values > 0);
    return OutputElementIterator(OutputElement(offsets_->n_values, shared_from_this()));
}

void OutputMeshBase::select_error_control_field(FieldSet &output_fields)
{
    if(refine_by_error_)
    {
        std::string error_control_field_name = "";
        // Read optional error control field name
        auto it = input_record_.find<std::string>("error_control_field");
        if(it) error_control_field_name = *it;

        FieldCommon* field =  output_fields.field(error_control_field_name);
        // throw input exception if the field is unknown
        if(field == nullptr){
            THROW(FieldSet::ExcUnknownField()
                    << FieldCommon::EI_Field(error_control_field_name)
                    << input_record_.ei_address());
            refine_by_error_ = false;
            return;
        }
        
        // throw input exception if the field is not scalar
        if( typeid(*field) == typeid(Field<3,FieldValue<3>::Scalar>) ) {
            
            error_control_field_ = static_cast<Field<3,FieldValue<3>::Scalar>*>(field);
            DebugOut() << "Output mesh will be refined according to field " << error_control_field_name << ".";
        }
        else{
            THROW(ExcFieldNotScalar()
                    << FieldCommon::EI_Field(error_control_field_name)
                    << input_record_.ei_address());
            refine_by_error_ = false;
            return;
        }
    }
    else
    {
        error_control_field_ = nullptr;
    }
}

bool OutputMeshBase::is_refined()
{
    return is_refined_;
}

unsigned int OutputMeshBase::n_elements()
{
    ASSERT_PTR(offsets_);
    return offsets_->n_values;
}

unsigned int OutputMeshBase::n_nodes()
{
    ASSERT_PTR(nodes_);
    return nodes_->n_values;
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


void OutputMesh::create_identical_mesh()
{
	DebugOut() << "Create outputmesh identical to computational one.";

    const unsigned int n_elements = orig_mesh_->n_elements(),
                       n_nodes = orig_mesh_->n_nodes();

    nodes_->data_.resize(spacedim*n_nodes);    // suppose 3D coordinates
    nodes_->n_values = n_nodes;
    unsigned int coord_id = 0,  // coordinate id in vector
                 node_id = 0;   // node id
    FOR_NODES(orig_mesh_, node) {
        node->aux = node_id;   // store node index in the auxiliary variable

        nodes_->data_[coord_id] = node->getX();  coord_id++;
        nodes_->data_[coord_id] = node->getY();  coord_id++;
        nodes_->data_[coord_id] = node->getZ();  coord_id++;
        node_id++;
    }
    
    orig_element_indices_ = std::make_shared<std::vector<unsigned int>>(n_elements);
    connectivity_->data_.reserve(4*n_elements);  //reserve - suppose all elements being tetrahedra (4 nodes)
    offsets_->data_.resize(n_elements);
    offsets_->n_values = n_elements;
    Node* node;
    unsigned int ele_id = 0,
                 connect_id = 0,
                 offset = 0,    // offset of node indices of element in node vector
                 li;            // local node index
    FOR_ELEMENTS(orig_mesh_, ele) {
        FOR_ELEMENT_NODES(ele, li) {
            node = ele->node[li];
            connectivity_->data_.push_back(node->aux);
            connect_id++;
        }
        
        // increase offset by number of nodes of the simplicial element
        offset += ele->dim() + 1;
        offsets_->data_[ele_id] = offset;
        (*orig_element_indices_)[ele_id] = ele_id;
        ele_id++;
    }
    connectivity_->data_.shrink_to_fit();
    connectivity_->n_values = connect_id;
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


void OutputMeshDiscontinuous::create_mesh(shared_ptr< OutputMesh > output_mesh)
{
    ASSERT_DBG(output_mesh->nodes_->n_values > 0);   //continuous data already computed
    
    if(nodes_->data_.size() > 0) return;          //already computed
    
    DebugOut() << "Create discontinuous outputmesh.";
    
    // connectivity = for every element list the nodes => its length corresponds to discontinuous data
    const unsigned int n_corners = output_mesh->connectivity_->n_values;

    // these are the same as in continuous case, so we copy the pointer.
    offsets_ = output_mesh->offsets_;
    orig_element_indices_ = output_mesh->orig_element_indices_;
    
    nodes_->data_.resize(spacedim*n_corners);
    nodes_->n_values = n_corners;
 
    connectivity_->data_.resize(n_corners);
    connectivity_->n_values = n_corners;

    unsigned int coord_id = 0,  // coordinate id in vector
                 corner_id = 0, // corner index (discontinous node)
                 li;            // local node index

    for(const auto & ele : *output_mesh)
    {
        unsigned int n = ele.n_nodes(), 
                     ele_idx = ele.idx(),
                     con_off = (* offsets_)[ele_idx];
                     
        for(li = 0; li < n; li++)
        {
            // offset of the first coordinate of the first node of the element in nodes_ vector
            unsigned int off = spacedim * (* output_mesh->connectivity_)[con_off - n + li];
            auto &d = output_mesh->nodes_->data_;
            
            nodes_->data_[coord_id] = d[off];   ++coord_id;
            nodes_->data_[coord_id] = d[off+1]; ++coord_id;
            nodes_->data_[coord_id] = d[off+2]; ++coord_id;
            
            connectivity_->data_[corner_id] = corner_id;
            corner_id++;
        }
    }
}


bool OutputMeshDiscontinuous::refinement_criterion()
{
    ASSERT(0).error("Not implemented yet.");
    return false;
}

void OutputMeshDiscontinuous::create_refined_mesh()
{
    nodes_ = std::make_shared<MeshData<double>>("", OutputDataBase::N_VECTOR);
    connectivity_ = std::make_shared<MeshData<unsigned int>>("connectivity");
    offsets_ = std::make_shared<MeshData<unsigned int>>("offsets");
    orig_element_indices_ = std::make_shared<std::vector<unsigned int>>();
    
    // index of last node added; set at the end of original ones
    unsigned int last_offset = 0;
    
    DBGMSG("start refinement\n");
    FOR_ELEMENTS(orig_mesh_, ele) {
        const unsigned int 
            dim = ele->dim(),
            ele_idx = ele->index();
        DBGMSG("ele index %d\n",ele_idx);
        
        AuxElement aux_ele;
        aux_ele.nodes.resize(ele->n_nodes());
        aux_ele.level = 0;
        
        Node* node; unsigned int li;
        FOR_ELEMENT_NODES(ele, li) {
            node = ele->node[li];
            aux_ele.nodes[li] = node->point();
        }
        
        std::vector<AuxElement> refinement;
        
        switch(dim){
            case 1: this->refine_aux_element<1>(aux_ele, refinement, ele->element_accessor(), error_control_field_); break;
            case 2: this->refine_aux_element<2>(aux_ele, refinement, ele->element_accessor(), error_control_field_); break;
            case 3: this->refine_aux_element<3>(aux_ele, refinement, ele->element_accessor(), error_control_field_); break;
            default: ASSERT(0 < dim && dim < 4);
        }
        
        //skip unrefined element
//         if(refinement.size() < 2) continue;
        
        unsigned int node_offset = nodes_->data_.size(),
                     con_offset = connectivity_->data_.size();
        nodes_->data_.resize(nodes_->data_.size() + (refinement.size() * (dim+1))*spacedim);
        connectivity_->data_.resize(connectivity_->data_.size() + refinement.size()*(dim+1));
//         orig_element_indices_->resize(orig_element_indices_->size() + refinement.size()*(dim+1));
        
        DBGMSG("ref size = %d\n", refinement.size());
        //gather coords and connectivity (in a continous way inside element)
        for(unsigned int i=0; i < refinement.size(); i++)
        {
            last_offset += dim+1;
            offsets_->data_.push_back(last_offset);
            (*orig_element_indices_).push_back(ele_idx);
            for(unsigned int j=0; j < dim+1; j++)
            {
                unsigned int con = i*(dim+1) + j;
                connectivity_->data_[con_offset + con] = con_offset + con;
                               
                for(unsigned int k=0; k < spacedim; k++) {
                    nodes_->data_[node_offset + con*spacedim + k] = refinement[i].nodes[j][k];
                }
            }
        }
    }
    
    connectivity_->n_values = connectivity_->data_.size();
    nodes_->n_values = nodes_->data_.size() / spacedim;
    offsets_->n_values = offsets_->data_.size();
    
    is_refined_ = true;
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
                                                 const ElementAccessor<spacedim> &ele_acc,
                                                 Field<3, FieldValue<3>::Scalar> *error_control_field )
{
    static const unsigned int n_subelements = 1 << dim;  //2^dim
    
    // FIXME Use RefElement<>::interact<> from intersections
    static const std::vector<std::vector<unsigned int>> line_nodes[] = {
        {}, //0D
        
        {{0,1}},
        
        {{0,1},
         {0,2},
         {1,2}},
        
        {{0,1},
         {0,2},
         {1,2},
         {0,3},
         {1,3},
         {2,3}}
    };
    
    static const std::vector<unsigned int> conn[] = {
        {}, //0D
        
        {0, 2,
         2, 1},
         
        {0, 3, 4,
         3, 1, 5,
         4, 5, 2,
         3, 5, 4},
         
        {0, 4, 5, 7,
         4, 1, 6, 8,
         5, 6, 2, 9,
         7, 8, 9, 3,
         4, 7, 8, 9,
         4, 6, 8, 9,
         4, 7, 5, 9,
         4, 5, 6, 9}
    }; 
//     DBGMSG("level = %d, %d\n", aux_element.level, max_refinement_level_);
 
    ASSERT_DBG(dim == aux_element.nodes.size()-1);
    
    // check refinement criterion
    bool is_refined_enough;
    if(refine_by_error_)
    {
        // compute centre of aux element
        Space<spacedim>::Point centre({0,0,0});
        for(auto& v : aux_element.nodes ) centre += v;
        centre = centre/aux_element.nodes.size();
        is_refined_enough = ! refinement_criterion_error(aux_element, centre, ele_acc, error_control_field);
    }   
    else
        is_refined_enough = ! refinement_criterion_uniform(aux_element);
    
    // if not refining any further, push into final vector
    if( is_refined_enough ) {
        refinement.push_back(aux_element);
        return;
    }
    
    std::vector<AuxElement> subelements(n_subelements);
    
    // FIXME Use RefElement<>::n_nodes and RefElement<>::n_lines from intersections
    const unsigned int n_old_nodes = dim+1,
                       n_new_nodes = (unsigned int)((dim * (dim + 1)) / 2); // new points are in the center of lines
    
    // auxiliary vectors
    std::vector<Space<spacedim>::Point> nodes = aux_element.nodes;
//     std::vector<unsigned int> node_numbering = aux_element.connectivity;
    nodes.reserve(n_old_nodes+n_new_nodes);
//     node_numbering.reserve(n_old_nodes+n_new_nodes);

    // create new points in the element
    for(unsigned int e=0; e < n_new_nodes; e++)
    {
        Space<spacedim>::Point p = nodes[line_nodes[dim][e][0]]+nodes[line_nodes[dim][e][1]];
        nodes.push_back( p / 2.0);
        //nodes.back().print();
        
//         last_node_idx++;
//         node_numbering.push_back(last_node_idx);
    }
   
    
    for(unsigned int i=0; i < n_subelements; i++)
    {
        AuxElement& sub_ele = subelements[i];
        sub_ele.nodes.resize(n_old_nodes);
//         sub_ele.connectivity.resize(n_old_nodes);
        sub_ele.level = aux_element.level+1;
        
        // over nodes
        for(unsigned int j=0; j < n_old_nodes; j++)
        {
            unsigned int conn_id = (n_old_nodes)*i + j;
            sub_ele.nodes[j] = nodes[conn[dim][conn_id]];
        }
        
        refine_aux_element<dim>(sub_ele, refinement, ele_acc, error_control_field);
    }
}

template void OutputMeshDiscontinuous::refine_aux_element<1>(const OutputMeshDiscontinuous::AuxElement&,std::vector< OutputMeshDiscontinuous::AuxElement >&, const ElementAccessor<spacedim> &, Field<3, FieldValue<3>::Scalar> *);
template void OutputMeshDiscontinuous::refine_aux_element<2>(const OutputMeshDiscontinuous::AuxElement&,std::vector< OutputMeshDiscontinuous::AuxElement >&, const ElementAccessor<spacedim> &, Field<3, FieldValue<3>::Scalar> *);
template void OutputMeshDiscontinuous::refine_aux_element<3>(const OutputMeshDiscontinuous::AuxElement&,std::vector< OutputMeshDiscontinuous::AuxElement >&, const ElementAccessor<spacedim> &, Field<3, FieldValue<3>::Scalar> *);

bool OutputMeshDiscontinuous::refinement_criterion_uniform(const OutputMeshDiscontinuous::AuxElement& ele)
{
    return (ele.level < max_level_);
}

bool OutputMeshDiscontinuous::refinement_criterion_error(const OutputMeshDiscontinuous::AuxElement& ele,
                                            const Space<spacedim>::Point &centre,
                                            const ElementAccessor<spacedim> &ele_acc,
                                            Field<3, FieldValue<3>::Scalar> *error_control_field
                                           )
{
    if(ele.level  < max_level_)
    {
        std::vector<double> nodes_val(ele.nodes.size());
        error_control_field->value_list(ele.nodes, ele_acc, nodes_val);
        
        double average_val = 0.0;
        for(double& v: nodes_val) 
            average_val += v;
        average_val = average_val / ele.nodes.size();
        
        double centre_val = error_control_field->value(centre,ele_acc);
        double diff = std::abs((average_val - centre_val)/centre_val);
        DBGMSG("diff: %f  %f  %f\n", diff, average_val, centre_val);
        return ( diff > refinement_error_tolerance_);
    }
    else
        return false;
}
