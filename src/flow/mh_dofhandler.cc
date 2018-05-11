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
 * @file    mh_dofhandler.cc
 * @brief   
 */

#include "flow/mh_dofhandler.hh"
#include "la/local_to_global_map.hh"
#include "mesh/long_idx.hh"
#include "mesh/mesh.h"
#include "mesh/partitioning.hh"
#include "mesh/side_impl.hh"
#include "system/sys_profiler.hh"

#include "intersection/mixed_mesh_intersections.hh"
#include "intersection/intersection_local.hh"

#include "fem/singularity.hh"
#include "fem/xfem_element_data.hh"
#include <fem/mapping_p1.hh>

const int MH_DofHandler::empty_node_idx = -1;

MH_DofHandler::MH_DofHandler()
:
el_4_loc(nullptr),
row_4_el(nullptr),
side_id_4_loc(nullptr),
side_row_4_id(nullptr),
edge_4_loc(nullptr),
row_4_edge(nullptr),
edge_ds(nullptr),
el_ds(nullptr),
side_ds(nullptr),
row_4_sing(nullptr),
row_4_vel_sing(nullptr),
row_4_press_sing(nullptr),
enrich_velocity(false),
enrich_pressure(false),
continuous_pu(false),
single_enr(false)
{}

MH_DofHandler::~MH_DofHandler()
{
    delete edge_ds;
    //  delete el_ds;
    delete side_ds;

    //  xfree(el_4_loc);
    delete [] row_4_el;
    delete [] side_id_4_loc;
    delete [] side_row_4_id;
    delete [] edge_4_loc;
    delete []  row_4_edge;
    delete [] row_4_sing;
    if(row_4_vel_sing) delete []  row_4_vel_sing;
    if(row_4_press_sing) delete []  row_4_press_sing;
}


void MH_DofHandler::reinit(Mesh *mesh) {
    mesh_ = mesh;
    elem_side_to_global.resize(mesh->n_elements() );
    FOR_ELEMENTS(mesh, ele) elem_side_to_global[ele.index()].resize(ele->n_sides());

    unsigned int i_side_global=0;
    FOR_ELEMENTS(mesh, ele) {
        for(unsigned int i_lside=0; i_lside < ele->n_sides(); i_lside++)
            elem_side_to_global[ele.index()][i_lside] = i_side_global++;
    }

    prepare_parallel();
    
    // convert row_4_id arrays from separate numberings to global numbering of rows
    make_row_numberings();
    
    // set dof offsets:
    offset_velocity = 0;
    offset_enr_velocity = mesh_->n_sides_;
    
    offset_pressure = offset_enr_velocity;
    offset_enr_pressure = offset_pressure + mesh_->n_elements();
    
    offset_edges = offset_enr_pressure;
    offset_enr_lagrange = offset_edges + mesh_->n_edges();
    
    print_dofs_dbg();
}


void MH_DofHandler::reinit(Mesh *mesh,
                           Field<3, FieldValue<3>::Scalar>& cross_section,
                           Field<3, FieldValue<3>::Scalar>& sigma) {
    ASSERT(enrich_pressure || enrich_velocity);
    
    mesh_ = mesh;
    elem_side_to_global.resize(mesh->n_elements() );
    FOR_ELEMENTS(mesh, ele) elem_side_to_global[ele.index()].resize(ele->n_sides());

    unsigned int i_side_global=0;
    FOR_ELEMENTS(mesh, ele) {
        for(unsigned int i_lside=0; i_lside < ele->n_sides(); i_lside++)
            elem_side_to_global[ele.index()][i_lside] = i_side_global++;
    }

    prepare_parallel();
    
//     // convert row_4_id arrays from separate numberings to global numbering of rows
//     make_row_numberings();
    if(enrich_velocity || enrich_pressure){
        if(xfem_dim == 2)
            create_enrichment(singularities_12d_, xfem_data_2d, cross_section, sigma);
        if(xfem_dim == 3)
            create_enrichment(singularities_13d_, xfem_data_3d, cross_section, sigma);
        
        // distribute FE enriched dofs
        distribute_enriched_dofs();
        // correct standard dofs:
        update_standard_dofs();
    }
    
    //HACK for a single processor
    unsigned int rows_starts = total_size();
    rows_ds = std::make_shared<Distribution>(&rows_starts, PETSC_COMM_WORLD);
    rows_ds->view(cout);
    
    print_dofs_dbg();
}


void MH_DofHandler::print_dofs_dbg()
{
    DBGCOUT(<< "\noffset_velocity " << offset_velocity << "\n"
            << "offset_enr_velocity " << offset_enr_velocity << "\n"
            << "offset_pressure " << offset_pressure << "\n"
            << "offset_enr_pressure " << offset_enr_pressure << "\n"
            << "offset_edges " << offset_edges << "\n"
            << "offset_enr_lagrange " << offset_enr_lagrange << "\n"
            << "total_size " << total_size() << "\n"
    );
    
    for (unsigned int i_loc = 0; i_loc < el_ds->lsize(); i_loc++) {
        auto ele_ac = accessor(i_loc);
        int dofs[200];
        int ndofs = ele_ac.get_dofs(dofs);
        
        DBGCOUT("### DOFS ele " << ele_ac.ele_global_idx() << "   ");
        cout << "[";
        if(ele_ac.is_enriched()) cout << ele_ac.xfem_data_pointer()->n_enrichments();
        else cout << "0";
        cout << "]  ";
        
        for(int i =0; i < ndofs; i++){
            cout << dofs[i] << " ";
        }
        cout << "\n";
    }
}



// ====================================================================================
// - compute optimal edge partitioning
// - compute appropriate partitioning of elements and sides
// - make arrays: *_id_4_loc and *_row_4_id to allow parallel assembly of the MH matrix
// ====================================================================================
void MH_DofHandler::prepare_parallel() {

    START_TIMER("prepare parallel");

    LongIdx *loc_part; // optimal (edge,el) partitioning (local chunk)
    LongIdx *id_4_old; // map from old idx to ids (edge,el)
    int loc_i;

    int e_idx;


    //ierr = MPI_Barrier(PETSC_COMM_WORLD);
    //OLD_ASSERT(ierr == 0, "Error in MPI_Barrier.");

    // row_4_el will be modified so we make a copy of the array from mesh
    row_4_el = new LongIdx[mesh_->n_elements()];
    std::copy(mesh_->get_row_4_el(), mesh_->get_row_4_el()+mesh_->n_elements(), row_4_el);
    el_4_loc = mesh_->get_el_4_loc();
    el_ds = mesh_->get_el_ds();

    //optimal element part; loc. els. id-> new el. numbering
    Distribution init_edge_ds(DistributionLocalized(), mesh_->n_edges(), PETSC_COMM_WORLD);
    // partitioning of edges, edge belongs to the proc of his first element
    // this is not optimal but simple
    loc_part = new LongIdx[init_edge_ds.lsize()];
    id_4_old = new LongIdx[mesh_->n_edges()];
    {
        loc_i = 0;
        FOR_EDGES(mesh_, edg ) {
            unsigned int i_edg = edg - mesh_->edges.begin();
            // partition
            e_idx = mesh_->element.index(edg->side(0)->element());
            if (init_edge_ds.is_local(i_edg)) {
                // find (new) proc of the first element of the edge
                loc_part[loc_i++] = el_ds->get_proc(row_4_el[e_idx]);
            }
            // id array
            id_4_old[i_edg] = i_edg;
        }
    }
    Partitioning::id_maps(mesh_->n_edges(), id_4_old, init_edge_ds, loc_part, edge_ds, edge_4_loc, row_4_edge);
    delete[] loc_part;
    delete[] id_4_old;

    // create map from mesh global edge id to new local edge id
    unsigned int loc_edge_idx=0;
    for (unsigned int i_el_loc = 0; i_el_loc < el_ds->lsize(); i_el_loc++) {
        auto ele = mesh_->element(el_4_loc[i_el_loc]);
        for (unsigned int i = 0; i < ele->n_sides(); i++) {
            unsigned int mesh_edge_idx= ele->side(i)->edge_idx();
            if ( edge_new_local_4_mesh_idx_.count(mesh_edge_idx) == 0 )
                // new local edge
                edge_new_local_4_mesh_idx_[mesh_edge_idx] = loc_edge_idx++;
        }
    }

    //optimal side part; loc. sides; id-> new side numbering
    Distribution init_side_ds(DistributionBlock(), mesh_->n_sides(), PETSC_COMM_WORLD);
    // partitioning of sides follows elements
    loc_part = new LongIdx[init_side_ds.lsize()];
    id_4_old = new LongIdx[mesh_->n_sides()];
    {
        int is = 0;
        loc_i = 0;
        FOR_SIDES(mesh_, side ) {
            // partition
            if (init_side_ds.is_local(is)) {
                // find (new) proc of the element of the side
                loc_part[loc_i++] = el_ds->get_proc(
                        row_4_el[mesh_->element.index(side->element())]);
            }
            // id array
            id_4_old[is++] = side_dof( side );
        }
    }

    Partitioning::id_maps(mesh_->n_sides(), id_4_old, init_side_ds, loc_part, side_ds,
            side_id_4_loc, side_row_4_id);
    delete [] loc_part;
    delete [] id_4_old;
}

// ========================================================================
// to finish row_4_id arrays we have to convert individual numberings of
// sides/els/edges to whole numbering of rows. To this end we count shifts
// for sides/els/edges on each proc and then we apply them on row_4_id
// arrays.
// we employ macros to avoid code redundancy
// =======================================================================
void MH_DofHandler::make_row_numberings() {
    int i, shift;
    int np = edge_ds->np();
    int edge_shift[np], el_shift[np], side_shift[np];
    unsigned int rows_starts[np];

    int edge_n_id = mesh_->n_edges(),
            el_n_id = mesh_->element.size(),
            side_n_id = mesh_->n_sides();

    // compute shifts on every proc
    shift = 0; // in this var. we count new starts of arrays chunks
    for (i = 0; i < np; i++) {
        side_shift[i] = shift - (side_ds->begin(i)); // subtract actual start of the chunk
        shift += side_ds->lsize(i);
        el_shift[i] = shift - (el_ds->begin(i));
        shift += el_ds->lsize(i);
        edge_shift[i] = shift - (edge_ds->begin(i));
        shift += edge_ds->lsize(i);
        rows_starts[i] = shift;
    }
    // apply shifts
    for (i = 0; i < side_n_id; i++) {
    	LongIdx &what = side_row_4_id[i];
        if (what >= 0)
            what += side_shift[side_ds->get_proc(what)];
    }
    for (i = 0; i < el_n_id; i++) {
    	LongIdx &what = row_4_el[i];
        if (what >= 0)
            what += el_shift[el_ds->get_proc(what)];

    }
    for (i = 0; i < edge_n_id; i++) {
    	LongIdx &what = row_4_edge[i];
        if (what >= 0)
            what += edge_shift[edge_ds->get_proc(what)];
    }
    // make distribution of rows
    for (i = np - 1; i > 0; i--)
        rows_starts[i] -= rows_starts[i - 1];

    rows_ds = std::make_shared<Distribution>(&(rows_starts[0]), PETSC_COMM_WORLD);
}


void MH_DofHandler::prepare_parallel_bddc() {
#ifdef FLOW123D_HAVE_BDDCML
    // auxiliary
    Element *el;
    LongIdx side_row, edge_row;

    global_row_4_sub_row = std::make_shared<LocalToGlobalMap>(rows_ds);

    //
    // ordering of dofs
    // for each subdomain:
    // | velocities (at sides) | pressures (at elements) | L. mult. (at edges) |
    for (unsigned int i_loc = 0; i_loc < el_ds->lsize(); i_loc++) {
        el = mesh_->element(el_4_loc[i_loc]);
        LongIdx el_row = row_4_el[el_4_loc[i_loc]];

        global_row_4_sub_row->insert( el_row );

        unsigned int nsides = el->n_sides();
        for (unsigned int i = 0; i < nsides; i++) {
            side_row = side_row_4_id[ side_dof( el->side(i) ) ];
            edge_row = row_4_edge[el->side(i)->edge_idx()];

            global_row_4_sub_row->insert( side_row );
            global_row_4_sub_row->insert( edge_row );
        }

        for (unsigned int i_neigh = 0; i_neigh < el->n_neighs_vb; i_neigh++) {
            // mark this edge
            edge_row = row_4_edge[el->neigh_vb[i_neigh]->edge_idx() ];
            global_row_4_sub_row->insert( edge_row );
        }
    }
    global_row_4_sub_row->finalize();
#endif // FLOW123D_HAVE_BDDCML
}



unsigned int MH_DofHandler::side_dof(const SideIter side) const {
    return elem_side_to_global[ side->element().index() ][ side->el_idx() ];
}


void MH_DofHandler::set_solution( double time, double * solution, double precision) {
	OLD_ASSERT( solution != NULL, "Empty solution.\n");
    mh_solution = solution;
    solution_precision = precision;
    time_ = time;
}

/// temporary replacement for DofHandler accessor, flux through given side
double MH_DofHandler::side_flux(const Side &side) const {
//     DBGVAR(elem_side_to_global[ side.element().index() ][ side.el_idx() ]);
    return mh_solution[ offset_velocity + elem_side_to_global[ side.element().index() ][ side.el_idx() ] ];
}

/// temporary replacement for DofHandler accessor, scalar (pressure) on edge of the side
double MH_DofHandler::side_scalar(const Side &side) const {
    unsigned int i_edg = side.edge_idx();
//     return mh_solution[ side.mesh()->n_sides() + side.mesh()->n_elements() + i_edg ];
    return mh_solution[ offset_edges + i_edg ];
}


double MH_DofHandler::element_scalar( ElementFullIter &ele ) const {
//     return mh_solution[ ele->mesh_->n_sides() + ele.index() ];
    return mh_solution[ offset_pressure + ele.index() ];
}


LocalElementAccessorBase<3> MH_DofHandler::accessor(uint local_ele_idx) {
    return LocalElementAccessorBase<3>(this, local_ele_idx);
}





/*************************************************************************************************************
 * ***********************************************************************************************************
 * **********************************************************************************************************/

void MH_DofHandler::clear_mesh_flags()
{
    mesh_flags_.resize(mesh_->n_elements());
    for(auto &&t : mesh_flags_)
        t = false;
}

void MH_DofHandler::clear_node_aux()
{
    FOR_ELEMENTS(mesh_, ele){
        for(unsigned int i=0; i < ele->n_nodes(); i++){
            ele->node[i]->aux = empty_node_idx;
        }
    }
}

unsigned int MH_DofHandler::n_enrichments()
{
    return (singularities_12d_.size() + singularities_13d_.size());
}

int MH_DofHandler::total_size()
{
    return offset_enr_lagrange + n_enrichments();
}

void MH_DofHandler::update_standard_dofs()
{
    
//     print_array(row_4_el, mesh_->n_elements(), "row_4_el(pressure)");
//     print_array(row_4_edge, mesh_->n_edges(), "row_4_edge(lagrange pressure)");
    
    unsigned int dof;
    dof = offset_pressure;
    for(unsigned int i=0; i < mesh_->n_elements(); i++, dof++){
        row_4_el[i] = dof;
    }
    
    dof = offset_edges;
    for(unsigned int i=0; i < mesh_->n_edges(); i++, dof++){
        row_4_edge[i] = dof;
    }
    
    dof = offset_enr_lagrange;
    
    unsigned int enr_size = n_enrichments();
    row_4_sing = new int[enr_size];
    for(unsigned int i=0; i < enr_size; i++, dof++)
        row_4_sing[i] = dof;

//     print_array(row_4_el, mesh_->n_elements(), "row_4_el(pressure)");
//     print_array(row_4_edge, mesh_->n_edges(), "row_4_edge(lagrange pressure)");
//     print_array(row_4_vel_sing, n_enrichments(), "row_4_vel_sing");
//     print_array(row_4_press_sing, n_enrichments(), "row_4_press_sing");
    
}


template<>
void MH_DofHandler::create_testing_singularities<Singularity<0>>(std::vector< Singularity0DPtr >& singularities,
                                                                 int & new_enrich_node_idx)
{
    //create singularity
//     Space<3>::Point center ({4.42, 2.08, 0}); //triangle
    Space<3>::Point center ({3.33, 3.33, 0}); //triangle
//     Space<3>::Point center ({-0.1, 0, 0}); //circle
//     Space<3>::Point center ({0.3, 0.3, 0}); //circle
    double radius = 0.03,//0.001643168,
           sigma_const = 10.0,
//            pressure = 50;//165.62;//100 + 28.5 * std::log(10);
           pressure = 100;
    Space<3>::Point direction_vector ({0,0,1});
    Space<3>::Point n;
    unsigned int n_qpoints = 1000;
    
    bool found = false;
    ElementFullIter ele2d(mesh_->element.begin());
    for(; ele2d != mesh_->element.end(); ++ele2d){
        if(ele2d->dim() == 2){
            MappingP1<2,3> map;
            arma::vec p = map.project_real_to_unit(center,map.element_map(*ele2d));
            if(map.is_point_inside(p))
            {
                n = arma::cross(ele2d->node[1]->point() - ele2d->node[0]->point(),
                                ele2d->node[2]->point() - ele2d->node[0]->point());
                found = true;
                break;
            }
        }
    }
    if(! found) return;
                
    DBGCOUT("singularity: 2d: " << ele2d->index() << "\n");
                
    auto sing = std::make_shared<Singularity<0>>(center, radius, direction_vector, n, n_qpoints);
    // set sigma of 1d element
    sing->set_sigma(sigma_const);
    sing->set_pressure(pressure);
    singularities.push_back(sing);

    //TODO: suggest proper enrichment radius
//     double enr_radius = 1.3*std::sqrt(ele2d->measure());
//     double enr_radius = 2.0;
    DBGCOUT(<< "enr_radius: " << enr_radius << "\n");
    clear_mesh_flags();

//     find_ele_to_enrich(singularities.back(), 0, ele2d, enr_radius, new_enrich_node_idx);
    enrich_ele(sing, singularities_12d_.size()-1, xfem_data_2d, 0, ele2d, new_enrich_node_idx);
    
    //flag the element
    mesh_flags_[ele2d->index()] = true;
    for(unsigned int n=0; n < ele2d->n_sides(); n++) {
        Edge* edge = ele2d->side(n)->edge();
        for(int j=0; j < edge->n_sides;j++) {
            if (edge->side(j)->element() != ele2d && edge->side(j)->element()->dim() == 2){
                find_ele_to_enrich(sing,0,edge->side(j)->element(),enr_radius, new_enrich_node_idx);
            }
        }
    }
    
    //create singularity
//     Space<3>::Point center2 ({2.08, 4.42, 0}); //triangle
//     double radius2 = 0.03,//0.001643168,
//            sigma_const2 = 10.0,
//            pressure2 = 50;//165.62;//100 + 28.5 * std::log(10);
//     Space<3>::Point direction_vector2 ({0,0,1});
//     Space<3>::Point n2;
//     
//     ele2d = mesh_->element.begin();
//     for(; ele2d != mesh_->element.end(); ++ele2d){
//         if(ele2d->dim() == 2){
//             MappingP1<2,3> map;
//             arma::vec p = map.project_real_to_unit(center2,map.element_map(*ele2d));
//             if(map.is_point_inside(p))
//             {
//                 n2 = arma::cross(ele2d->node[1]->point() - ele2d->node[0]->point(),
//                                 ele2d->node[2]->point() - ele2d->node[0]->point());
//                 break;
//             }
//         }
//     }
//                 
//     DBGCOUT("singularity: 2d: " << ele2d->index() << "\n");
//                 
//     auto sing2 = std::make_shared<Singularity0D>(center2, radius2, direction_vector2, n2);
//     // set sigma of 1d element
//     sing2->set_sigma(sigma_const2);
//     sing2->set_pressure(pressure2);
//     singularities.push_back(sing2);

//     clear_mesh_flags();
//     find_ele_to_enrich(singularities.back(), 0, ele2d, enr_radius, new_enrich_node_idx);
}


template<>
void MH_DofHandler::create_testing_singularities<Singularity<1>>(std::vector< Singularity1DPtr >& singularities,
                                                                int & new_enrich_node_idx)
{
    typedef Space<3>::Point Point;
    //create singularity
    Point a ({3.33, 3.33, 0});
    Point b ({3.33, 3.33, 10});
    double radius = 0.03,
           sigma_const = 10.0,
           pressure = 100;
    unsigned int n = 300, m = 1000;
    
    Point center = a;
    bool found = false;
    ElementFullIter ele3d(mesh_->element.begin());
    for(; ele3d != mesh_->element.end(); ++ele3d){
        if(ele3d->dim() == 3){
            MappingP1<3,3> map;
            arma::vec p = map.project_real_to_unit(center,map.element_map(*ele3d));
            if(map.is_point_inside(p)){
                found = true;
                break;
            }
        }
    }
    if(! found) return;
                
    DBGCOUT("singularity: 3d: " << ele3d->index() << "\n");
                
    auto sing = std::make_shared<Singularity<1>>(a,b,radius,n,m);
    // set sigma of 1d element
    sing->set_sigma(sigma_const);
    sing->set_pressure(pressure);
    singularities.push_back(sing);

    //TODO: suggest proper enrichment radius
    DBGCOUT(<< "enr_radius: " << enr_radius << "\n");
    clear_mesh_flags();

//     find_ele_to_enrich(singularities.back(), 0, ele3d, enr_radius, new_enrich_node_idx);
    enrich_ele(sing, singularities_13d_.size()-1, xfem_data_3d, 0, ele3d, new_enrich_node_idx);
    
    //flag the element
    mesh_flags_[ele3d->index()] = true;
    for(unsigned int n=0; n < ele3d->n_sides(); n++) {
        Edge* edge = ele3d->side(n)->edge();
        for(int j=0; j < edge->n_sides;j++) {
            if (edge->side(j)->element() != ele3d && edge->side(j)->element()->dim() == 3){
                find_ele_to_enrich(sing,0,edge->side(j)->element(),enr_radius, new_enrich_node_idx);
            }
        }
    }
}


template<>
std::shared_ptr<Singularity<0>> MH_DofHandler::create_sing<2>(IntersectionLocal<1,2>* il,
                                    double cross_section,
                                    double sigma)
{
    DebugOut().fmt("intersection:  c {} b {}\n", il->component_ele_idx(), il->bulk_ele_idx());
    ElementFullIter ele = mesh_->element(il->component_ele_idx());
    ElementFullIter bulk_ele = mesh_->element(il->bulk_ele_idx());
    //create singularity
    Space<3>::Point center = (*il)[0].coords(ele);
//     center.print(DebugOut(),"center");
    
    double radius = std::sqrt(cross_section/M_PI);
    DebugOut() << "radius " << radius << "\n";
    
    const unsigned int n_qpoints = 1000;
    Space<3>::Point n = arma::cross(bulk_ele->node[1]->point() - bulk_ele->node[0]->point(),
                                    bulk_ele->node[2]->point() - bulk_ele->node[0]->point());
    Space<3>::Point direction_vector(ele->node[1]->point() - ele->node[0]->point());
    
    auto sing = std::make_shared<Singularity<0>>(center, radius, direction_vector, n, n_qpoints);
    // set sigma of 1d element
    sing->set_sigma(sigma);
    return sing;
}

template<>
std::shared_ptr<Singularity<1>> MH_DofHandler::create_sing<3>(IntersectionLocal<1,3>* il,
                                    double cross_section,
                                    double sigma)
{
    DebugOut().fmt("intersection:  c {} b {}\n", il->component_ele_idx(), il->bulk_ele_idx());
    ElementFullIter ele = mesh_->element(il->component_ele_idx());
    
    //create singularity
    double radius = std::sqrt(cross_section/M_PI);
    DebugOut() << "radius " << radius << "\n";
    
    Space<3>::Point a = ele->node[0]->point();
    Space<3>::Point b = ele->node[1]->point();
//     a.print(DebugOut(),"pointA");
//     b.print(DebugOut(),"pointB");
    
    const unsigned int n_qpoints = 300;
    const unsigned int m_qpoints = 1000;
    
    auto sing = std::make_shared<Singularity<1>>(a,b, radius, n_qpoints, m_qpoints);
    // set sigma of 1d element
    sing->set_sigma(sigma);
    return sing;
}

template<int dim>
void MH_DofHandler::create_enrichment(std::vector<std::shared_ptr<Singularity<dim-2>>>& singularities,
                                      std::vector<XFEMElementSingularData<dim>>& xfem_data,
                                      Field<3, FieldValue<3>::Scalar>& cross_section,
                                      Field<3, FieldValue<3>::Scalar>& sigma)
{
    DBGCOUT(<< "MH_DofHandler - create_singularities " << dim << "d\n");
    //TODO:
    //- count different types of intersections inside InspectElements to be able to allocate other objects
    //- differ 1d-2d intersections: one point intersection in plane versus in 3D
    
    singularities.clear();
    clear_node_aux();
    
//     unsigned int sing_count_guess = 100;
//     singularities.reserve(sing_count_guess);
//     node_values.reserve(sing_count_guess);
//     node_vec_values.reserve(sing_count_guess);
    
    //TODO propose some allocation size for xfem data
    xfem_data.reserve(mesh_->n_elements());
    
    int new_enrich_node_idx = 0;
//     create_testing_singularities(singularities, new_enrich_node_idx);

    for (unsigned int i_loc = 0; i_loc < el_ds->lsize(); i_loc++) {
        auto ele_ac = accessor(i_loc);
        ElementFullIter ele = ele_ac.full_iter();
        if(ele_ac.dim() == 1) {
            auto &isec_list = mesh_->mixed_intersections().element_intersections_[ele_ac.ele_global_idx()];
            if(isec_list.size() == 0) continue;
            
            for(auto &isec : isec_list ) {
                ASSERT_PTR_DBG(isec.second);
                ElementFullIter bulk_ele = mesh_->element(isec.second->bulk_ele_idx());
                
                Space<3>::Point center = ele->centre();
                IntersectionLocal<1,dim>* il = static_cast<IntersectionLocal<1,dim>*>(isec.second);
                
                if(bulk_ele->dim() == 2){
                    if(il->size() != 1) continue;   //process only IL with one IP
                    
                    center = (*il)[0].coords(ele);
                }
                
                DBGCOUT("singularity: ele " << dim << "d: " << ele->index() << "  bulk_ele: " << bulk_ele->index() << "\n");
                double cs = cross_section.value(center, ele->element_accessor()); // pi*r^2
                double sgm = sigma.value(center, ele->element_accessor());
                
                singularities.push_back(create_sing<dim>(il, cs, sgm));
                
                DBGCOUT(<< "enr_radius: " << enr_radius << "\n");
                clear_mesh_flags();
                find_ele_to_enrich(singularities.back(), ele->index(), bulk_ele, enr_radius, new_enrich_node_idx);
                
//                 // in 3d, create only one singularity per 1d element
//                 if(dim == 3) break;
                //TODO: in 2d detect whether 2d intersection elements (of one 1d element) are lying in a common plane
                // if so, then create only single singularity, otherwise create another singularity in the other plane
                // => create a single singularity per detected plane
                // both in 2d and 3d, create only one singularity per 1d element
                break;
            }
        }
    }
    
//     FOR_ELEMENTS(mesh_,ele){
//         unsigned int idx = ele->index();
//         if(ele->dim() == 1) {
//             
//             std::vector<computeintersection::ILpair>& ilpairs = intersections->intersection_map_[idx];
//             if(ilpairs.size() == 0) continue;
//             
//             for(computeintersection::ILpair& ilp : ilpairs){
//                 ASSERT_PTR_DBG(ilp.second);
//                 //TODO: check different intersection
//                 computeintersection::IntersectionLocal<1,2>* il = static_cast<computeintersection::IntersectionLocal<1,2>*>(ilp.second);
//                 ASSERT_PTR_DBG(il);
//                 // we want to consider here only singularities (line segment and triangle not in a plane)
//                 if(il->size() != 1) continue;   //process only IL with one IP
//                 
//                 ElementFullIter ele2d = mesh_->element(il->bulk_ele_idx());
//                 
//                 DBGCOUT("singularity: ele 1d: " << idx << "  2d: " << ele2d->index() << "\n");
//                 //create singularity
//                 Space<3>::Point center = (*il)[0].coords(ele);
//                 double radius = cross_section.value(center,ele->element_accessor());
// //                 DebugOut().fmt("intersection:  c {} b {}\n", il->component_idx(), il->bulk_ele_idx());
// //                 center.print(MessageOut(),"center");
// //                 DebugOut() << "radius " << radius << "\n";
//                 Space<3>::Point n = arma::cross(ele2d->node[1]->point() - ele2d->node[0]->point(),
//                                                 ele2d->node[2]->point() - ele2d->node[0]->point());
//                 Space<3>::Point direction_vector(ele->node[1]->point() - ele->node[0]->point());
//                 
//                 auto sing = std::make_shared<Singularity0D<3>>(center, radius, direction_vector, n);
//                 // set sigma of 1d element
//                 sing->set_sigma(sigma.value(center, ele->element_accessor()));
//                 singularities.push_back(sing);
// //                 node_values.push_back(std::map<int, double>());
// //                 node_vec_values.push_back(std::map<int, Space<3>::Point>());
//                 
// //                 unsigned int sing_idx = singularities_12d_.size()-1;
// //                 if(ele->xfem_data == nullptr){
// //                     xfem_data_1d.push_back(XFEMComplementData(sing, sing_idx));
// //                     xfem_data_1d.back().set_element(idx, ele2d->index());
// //                     xfem_data_1d.back().set_complement();
// //                     ele->xfem_data = & xfem_data_1d.back();
// //                 }
// //                 else{
// //                     auto xdata = static_cast<XFEMComplementData*>(ele->xfem_data);
// //                     ASSERT_DBG(xdata != nullptr).error("XFEM data object is not of XFEMComplementData Type!");
// //                     xdata->add_data(sing, sing_idx);
// //                 }
//                 
//                 
//                 //TODO: suggest proper enrichment radius
//                 double enr_radius = 1.5*std::sqrt(ele2d->measure());
//                 DBGCOUT(<< "enr_radius: " << enr_radius << "\n");
//                 clear_mesh_flags();
// 
//                 find_ele_to_enrich(singularities.back(), idx, ele2d, enr_radius, new_enrich_node_idx);
//             }
//         }
//     }
    
    //shrink here (invalidates iterators, pointers); element xdata pointer is set in distribute_enriched_dofs()
    xfem_data.shrink_to_fit();

    // Evaluate singularity quad points.
    // Update invalid pointers.
//     if(! (enrich_pressure || enrich_velocity))
    for(auto& xdata : xfem_data){
        ElementFullIter ele = mesh_->element(xdata.ele_global_idx());
        xdata.create_sing_quads(ele);
//         DBGVAR(xdata.n_enrichments_intersect());
//         xdata.print(cout);
        ele->xfem_data = &xdata;
    }
    
//     singularities.shrink_to_fit();
//     node_values.shrink_to_fit();
//     node_vec_values.shrink_to_fit();
//     clear_mesh_flags();
}

void MH_DofHandler::distribute_enriched_dofs()
{
    // set dof offsets:
    offset_velocity = 0;
    offset_enr_velocity = mesh_->n_sides_;
    
    int temp_offset;
//     unsigned int max_enr_per_node = 1;
    
    unsigned int enr_size = singularities_12d_.size() + singularities_13d_.size();
    
    //distribute enriched dofs:
    temp_offset = offset_enr_velocity; // will return last dof + 1 (it means new available dof)
    if(enrich_velocity) {
        if(single_enr){
            row_4_vel_sing = new int[enr_size];
            for(unsigned int i=0; i < singularities_12d_.size(); i++, temp_offset++)
                row_4_vel_sing[i] = temp_offset;
            for(unsigned int i=singularities_12d_.size(); i < enr_size; i++, temp_offset++)
                row_4_vel_sing[i] = temp_offset;
        }
        else
        {
//             if(continuous_pu){
//                 // temporary dof vector
//                 std::vector<std::vector<int>> enr_dofs_velocity(n_enriched_nodes, std::vector<int>(max_enr_per_node, empty_node_idx));
//                 distribute_enriched_dofs(enr_dofs_velocity, temp_offset, Quantity::velocity);
//             }
//             else distribute_enriched_dofs(temp_offset, Quantity::velocity);
            distribute_enriched_dofs(temp_offset, Quantity::velocity);
        }
    }
    
    offset_pressure = temp_offset;
    offset_enr_pressure = offset_pressure + mesh_->n_elements();
    
    temp_offset = offset_enr_pressure; // will return last dof + 1 (it means new available dof)
    if(enrich_pressure){
        if(single_enr){
            row_4_press_sing = new int[enr_size];
            for(unsigned int i=0; i < singularities_12d_.size(); i++, temp_offset++)
                row_4_press_sing[i] = temp_offset;
            for(unsigned int i=singularities_12d_.size(); i < enr_size; i++, temp_offset++)
                row_4_press_sing[i] = temp_offset;
        }
        else{
//             if(continuous_pu){
                // temporary dof vector
//                 std::vector<std::vector<int>> enr_dofs_pressure(n_enriched_nodes, std::vector<int>(max_enr_per_node, empty_node_idx));
//                 distribute_enriched_dofs(enr_dofs_pressure, temp_offset, Quantity::pressure);
//             }
//             else distribute_enriched_dofs(temp_offset, Quantity::pressure);
            distribute_enriched_dofs(temp_offset, Quantity::pressure);
        }
    }
    offset_edges = temp_offset;
    offset_enr_lagrange = offset_edges + mesh_->n_edges();
//     offset_enr_lagrange = offset_edges + mesh_->n_edges() - 1;
    
}

// void MH_DofHandler::distribute_enriched_dofs(int n_enriched_nodes)
// {
//     // set dof offsets:
//     offset_velocity = 0;
//     offset_enr_velocity = mesh_->n_sides_;
//     
//     int temp_offset;
//     unsigned int max_enr_per_node = 1;
//     
//     //distribute enriched dofs:
//     temp_offset = offset_enr_velocity; // will return last dof + 1 (it means new available dof)
//     if(enrich_velocity) {
//         if(single_enr){
//             row_4_vel_sing = new int[singularities_12d_.size()];
//             for(unsigned int i=0; i < singularities_12d_.size(); i++, temp_offset++)
//                 row_4_vel_sing[i] = temp_offset;
//         }
//         else
//         {
//             if(continuous_pu){
//                 // temporary dof vector
//                 std::vector<std::vector<int>> enr_dofs_velocity(n_enriched_nodes, std::vector<int>(max_enr_per_node, empty_node_idx));
//                 distribute_enriched_dofs(enr_dofs_velocity, temp_offset, Quantity::velocity);
//             }
//             else distribute_enriched_dofs(temp_offset, Quantity::velocity);
//         }
//     }
//     
//     offset_pressure = temp_offset;
//     offset_enr_pressure = offset_pressure + mesh_->n_elements();
//     
//     temp_offset = offset_enr_pressure; // will return last dof + 1 (it means new available dof)
//     if(enrich_pressure){
//         if(single_enr){
//             row_4_press_sing = new int[singularities_12d_.size()];
//             for(unsigned int i=0; i < singularities_12d_.size(); i++, temp_offset++)
//                 row_4_press_sing[i] = temp_offset;
//         }
//         else{
//             if(continuous_pu){
//                 // temporary dof vector
//                 std::vector<std::vector<int>> enr_dofs_pressure(n_enriched_nodes, std::vector<int>(max_enr_per_node, empty_node_idx));
//                 distribute_enriched_dofs(enr_dofs_pressure, temp_offset, Quantity::pressure);
//             }
//             else distribute_enriched_dofs(temp_offset, Quantity::pressure);
//         }
//     }
//     offset_edges = temp_offset;
// //     offset_enr_lagrange = offset_edges + mesh_->n_edges();
//     offset_enr_lagrange = offset_edges + mesh_->n_edges() - 1;
//     
// }

template<int dim, class Enr>
void MH_DofHandler::enrich_ele(std::shared_ptr<Enr> sing,
                               unsigned int sing_idx,
                               std::vector<XFEMElementSingularData<dim>>& xfem_data,
                               int ele1d_global_idx,
                               ElementFullIter ele,
                               int& new_enrich_node_idx)
{
    DBGVAR(ele->index());
    
    XFEMElementSingularData<dim> * xdata;
    
    if(ele->xfem_data == nullptr){   //possibly create new one
        xfem_data.push_back(XFEMElementSingularData<dim>());
        xdata = & xfem_data.back();
//             xdata->set_node_values(&node_values, &node_vec_values);
        //TODO: set number of quantities
        xdata->global_enriched_dofs().resize(2);
        
        //HACK for xfem without enriching:
        if(! (enrich_pressure || enrich_velocity)){
            xdata->global_enriched_dofs()[0].resize(1);
            xdata->global_enriched_dofs()[1].resize(1);
        }
        xdata->set_element(ele->index());
        ele->xfem_data = xdata;
    }
    else{
//         DBGCOUT(<< "existing XData\n");
        xdata = static_cast<XFEMElementSingularData<dim>*>(ele->xfem_data);
//         xdata->print(cout);
        ASSERT_DBG(xdata != nullptr).error("XFEM data object is not of XFEMElementSingularData<dim> Type!");
    }
    
    xdata->add_data(sing, sing_idx, ele1d_global_idx);
    
    //HACK for xfem without enriching:
    // shortcut when not enriching
//         if(! (enrich_velocity || enrich_pressure)) return;
    
    Node* node; //shortcut
    // number the enriched nodes and compute node values
    for(unsigned int i=0; i < ele->n_nodes(); i++){
        node = ele->node[i];
        // number enriched nodes
        if (node->aux == empty_node_idx){
//                 DBGCOUT(<< "node number: " << new_enrich_node_idx << "\n");
            node->aux = new_enrich_node_idx;
            new_enrich_node_idx++;
        }
        
        // map.insert does not do anything if key already exists
//             node_values[sing_idx].insert(std::make_pair(node->aux, sing->value(node->point())) );       //pressure
//             node_vec_values[sing_idx].insert(std::make_pair(node->aux, sing->grad(node->point())) );    //velocity
    }
}


void MH_DofHandler::find_ele_to_enrich(Singularity0DPtr sing,
                                       int ele1d_global_idx,
                                   ElementFullIter ele,
                                   double radius,
                                   int& new_enrich_node_idx
                                  )
{   
//     DBGVAR(ele->index());
    // check flag at the element so element is checked only once
    if(mesh_flags_[ele->index()]) return;
    
    //flag the element
    mesh_flags_[ele->index()] = true;
//     bool enrich = true;
    bool enrich = false;
    for(unsigned int i=0; i < ele->n_nodes(); i++){
        double d = sing->geometry().distance(ele->node[i]->point());
//         DBGCOUT(<< d << "\n");
        if(d < radius){
            enrich = true;
        }
        else
        {
            MappingP1<2,3> map;
            arma::vec up = map.project_real_to_unit(sing->geometry_ellipse().center(),map.element_map(*ele));
            if(map.is_point_inside(up))
                enrich = true;
        }
    }
    
    // front advancing enrichment of neighboring elements
    if(enrich){
        enrich_ele(sing, singularities_12d_.size()-1, xfem_data_2d, ele1d_global_idx, ele, new_enrich_node_idx);
        
//         DebugOut() << "n_neighs_vb " << ele->n_neighs_vb << "\n";
        for(unsigned int n=0; n < ele->n_sides(); n++) {
            Edge* edge = ele->side(n)->edge();
            for(int j=0; j < edge->n_sides;j++) {
                if (edge->side(j)->element() != ele && edge->side(j)->element()->dim() == 2){
//                     DebugOut() << "Go to ele " << edge->side(j)->element()->index() << "\n";
                    find_ele_to_enrich(sing,ele1d_global_idx,edge->side(j)->element(),radius, new_enrich_node_idx);
                }
            }
        }
    }
}

void MH_DofHandler::find_ele_to_enrich(Singularity1DPtr sing,
                                       int ele1d_global_idx,
                                       ElementFullIter ele,
                                       double radius,
                                       int& new_enrich_node_idx)
{   
    typedef Space<3>::Point Point;
//     DBGVAR(ele->index());
    // check flag at the element so element is checked only once
    if(mesh_flags_[ele->index()]) return;
    
    //flag the element
    mesh_flags_[ele->index()] = true;
//     bool enrich = true;
    
    bool enrich = false;
    // enrich all other 3d elements intersecting the 1d
    auto& intersections = mesh_->mixed_intersections().element_intersections_[ele1d_global_idx];
    for(auto& il:intersections){
//         DBGCOUT(<< il.first << "\n");
        if(il.first == ele->index()){
            enrich = true;
            break;
        }
    }
    
    if(! enrich){
        const CylinderGeometry& geom = sing->geometry_cylinder();
        for(unsigned int i=0; i < ele->n_nodes(); i++){
            Point dp = geom.dist_vector(ele->node[i]->point());
            double d = arma::norm(dp,2);
            
            // test distance to line
            if(d < radius){
                
//                 DebugOut().fmt("d: {} < enr_r: {}\n",d,radius);
                Point p = ele->node[i]->point() - dp;
                MappingP1<1,3> map;
                arma::vec up = map.project_real_to_unit(p,map.element_map(mesh_->element[ele1d_global_idx]));
                if(map.is_point_inside(up)){
                    enrich = true;
                    break;
                }
                
                //test distance to A
                double da = arma::norm(geom.a()-ele->node[i]->point(),2);
                if(da < radius){
                    enrich = true;
                    break;
                }
                
                //test distance to B
                double db = arma::norm(geom.b()-ele->node[i]->point(),2);
                if(db < radius){
                    enrich = true;
                    break;
                }
            }
        }
    }
    
    // front advancing enrichment of neighboring elements
    if(enrich){
//         DebugOut() << "Enrich when intersecting.\n";
        enrich_ele(sing, singularities_13d_.size()-1, xfem_data_3d, ele1d_global_idx, ele, new_enrich_node_idx);
        
//         DebugOut() << "n_neighs_vb " << ele->n_neighs_vb << "\n";
        for(unsigned int n=0; n < ele->n_sides(); n++) {
            Edge* edge = ele->side(n)->edge();
            for(int j=0; j < edge->n_sides;j++) {
                if (edge->side(j)->element() != ele && edge->side(j)->element()->dim() == 3){
//                     DebugOut() << "Go to ele " << edge->side(j)->element()->index() << "\n";
                    find_ele_to_enrich(sing,ele1d_global_idx,edge->side(j)->element(),radius, new_enrich_node_idx);
                }
            }
        }
    }
}

void MH_DofHandler::distribute_enriched_dofs(vector< std::vector< int > >& temp_dofs,
                                             int& offset,
                                             Quantity quant)
{
    unsigned int w,i,node_idx;
    
    for(XFEMElementDataBase& xdata : xfem_data_2d){
        ElementFullIter ele = mesh_->element(xdata.ele_global_idx());
        std::vector<std::vector<int>>& dofs = xdata.global_enriched_dofs()[quant];
        dofs.resize(xdata.n_enrichments(), std::vector<int>(ele->n_nodes(), -1));
//         DBGCOUT(<<"dofs xdata\n");
        for(w=0; w < xdata.n_enrichments(); w++){
            for(i=0; i < ele->n_nodes(); i++){
                node_idx = ele->node[i]->aux;
                if(temp_dofs[node_idx][w] == empty_node_idx){
//                     DBGCOUT(<< node_idx << " new dof " << offset << "\n");
                    dofs[w][i] = temp_dofs[node_idx][w] = offset;
                    offset++;
                }
                else{
//                     DBGCOUT(<< node_idx << " old dof " << temp_dofs[node_idx][w] << "\n");
                    dofs[w][i] = temp_dofs[node_idx][w];
                }
            }
        }
//         xdata.print(cout);
        ele->xfem_data = &xdata;
        
//         DBGCOUT(<< "xd[0]: " <<  xdata.global_enriched_dofs()[0].size() << "\n");
    }
}

// void MH_DofHandler::distribute_enriched_dofs(int& offset,
//                                              Quantity quant)
// {
//     unsigned int w,i;
//     
//     for(XFEMElementSingularData& xdata : xfem_data){
//         ElementFullIter ele = mesh_->element(xdata.ele_global_idx());
//         std::vector<std::vector<int>>& dofs = xdata.global_enriched_dofs()[quant];
//         dofs.resize(xdata.n_enrichments(), std::vector<int>(ele->n_nodes(), -1));
// //         DBGCOUT(<<"dofs xdata\n");
//         for(w=0; w < xdata.n_enrichments(); w++){
//             for(i=0; i < ele->n_nodes(); i++){
//                 dofs[w][i] = offset;
//                 offset++;
//             }
//         }
// //         xdata.print(cout);
//         ele->xfem_data = &xdata;
//         
// //         DBGCOUT(<< "xd[0]: " <<  xdata.global_enriched_dofs()[0].size() << "\n");
//     }
// }


void MH_DofHandler::distribute_enriched_dofs(int& offset,
                                             Quantity quant)
{
    unsigned int w;
    
    for(XFEMElementDataBase& xdata : xfem_data_2d){
        ElementFullIter ele = mesh_->element(xdata.ele_global_idx());
        std::vector<std::vector<int>>& dofs = xdata.global_enriched_dofs()[quant];
        dofs.resize(xdata.n_enrichments(), std::vector<int>(1, -1));
//         DBGCOUT(<<"dofs xdata\n");
        for(w=0; w < xdata.n_enrichments(); w++){
            dofs[w][0] = offset;
            offset++;
        }
//         xdata.print(cout);
        ele->xfem_data = &xdata;
        
//         DBGCOUT(<< "xd[0]: " <<  xdata.global_enriched_dofs()[0].size() << "\n");
    }
}
