/*
 * inspectelements.cpp
 *
 *  Created on: 29.4.2016
 *      Author: pe
 */

#include "inspect_elements_algorithm.hh"
#include "intersection_point_aux.hh"
#include "intersection_aux.hh"
#include "intersection_local.hh"
#include "compute_intersection.hh"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "mesh/ref_element.hh"
#include "mesh/bih_tree.hh"

namespace computeintersection {

template<unsigned int simplex_dim>
void update_simplex(const ElementFullIter& element, Simplex< simplex_dim >& simplex)
{
    arma::vec3 *field_of_points[simplex_dim+1];
    for(unsigned int i=0; i < simplex_dim+1; i++)
        field_of_points[i]= &(element->node[i]->point());
    simplex.set_simplices(field_of_points);
}


template<unsigned int dim>    
InspectElementsAlgorithm<dim>::InspectElementsAlgorithm(Mesh* _mesh)
: mesh(_mesh)
{
}

template<unsigned int dim>   
InspectElementsAlgorithm<dim>::~InspectElementsAlgorithm()
{}


    
template<unsigned int dim>
void InspectElementsAlgorithm<dim>::init()
{
    START_TIMER("Intersection initialization");
    last_slave_for_3D_elements.assign(mesh->n_elements(), undefined_elm_idx_);
    closed_elements.assign(mesh->n_elements(), false);
    intersection_list_.assign(mesh->n_elements(),std::vector<IntersectionAux<dim,3>>());
    n_intersections_ = 0;
    component_counter_ = 0;
    END_TIMER("Intersection initialization");
}

template<unsigned int dim>
void InspectElementsAlgorithm<dim>::compute_bounding_boxes()
{
    START_TIMER("Compute bounding boxes");
    if(elements_bb.size() == 0){
        elements_bb.resize(mesh->n_elements());
        bool first_3d_element = true;
        FOR_ELEMENTS(mesh, elm) {

            elements_bb[elm->index()] = elm->bounding_box();

                if (elm->dim() == 3){
                    if(first_3d_element){
                        first_3d_element = false;
                        mesh_3D_bb = elements_bb[elm->index()];
                    }else{
                        mesh_3D_bb.expand(elements_bb[elm->index()].min());
                        mesh_3D_bb.expand(elements_bb[elm->index()].max());
                    }
                }
        }
    }
    END_TIMER("Compute bounding boxes");
}

template<unsigned int dim> 
bool InspectElementsAlgorithm<dim>::compute_initial_CI(unsigned int component_ele_idx,
                                                       unsigned int bulk_ele_idx,
                                                       unsigned int component_idx,
                                                       std::vector< unsigned int >& prolongation_table)
{
    IntersectionAux<dim,3> is(component_ele_idx, bulk_ele_idx, component_idx);
    START_TIMER("Compute intersection");
    ComputeIntersection<Simplex<dim>, Simplex<3>> CI(component_simplex, tetrahedron);
    CI.init();
    CI.compute(is, prolongation_table);
    END_TIMER("Compute intersection");
    
    last_slave_for_3D_elements[bulk_ele_idx] = component_ele_idx;
    
    if(is.points().size() > 0) {
        
        intersection_list_[component_ele_idx].push_back(is);
        n_intersections_++;
        return true;
    }
    else return false;
}




template<unsigned int dim>
bool InspectElementsAlgorithm<dim>::intersection_exists(unsigned int component_ele_idx, unsigned int bulk_ele_idx) 
{
    START_TIMER("intersection_exists");
    bool found = false;
    
    START_TIMER("cycle");
    unsigned int i = 0;
    for(; i < intersection_list_[component_ele_idx].size();i++){

        if(intersection_list_[component_ele_idx][i].bulk_ele_idx() == bulk_ele_idx){
            found = true;
            break;
        }
    }
    ADD_CALLS(i);
    END_TIMER("cycle");
    END_TIMER("intersection_exists");
    return found;
}


template<unsigned int dim>
void InspectElementsAlgorithm<dim>::compute_intersections(const BIHTree& bih)
{
    DebugOut() << "#########   ALGORITHM: compute_intersections   #########\n";
    
    init();
    
    START_TIMER("Element iteration");
    
    FOR_ELEMENTS(mesh, elm) {
        unsigned int component_ele_idx = elm->index();
        
        if (elm->dim() == dim &&                                // is component element
            !closed_elements[component_ele_idx] &&                    // is not closed yet
            bih.ele_bounding_box(component_ele_idx).intersect(bih.tree_box()))    // its bounding box intersects 3D mesh bounding box
        {    
            std::vector<unsigned int> searchedElements;
            
            START_TIMER("BIHtree find");
            bih.find_bounding_box(bih.ele_bounding_box(component_ele_idx), searchedElements);
            END_TIMER("BIHtree find");

//             component_counter_++;
//             DebugOut().fmt("comp: {}\n", component_counter_);
            
            START_TIMER("Bounding box element iteration");
            
            // Go through all element which bounding box intersects the component element bounding box
            for (std::vector<unsigned int>::iterator it = searchedElements.begin(); it!=searchedElements.end(); it++)
            {
                unsigned int bulk_ele_idx = *it;
                ElementFullIter ele_3D = mesh->element(bulk_ele_idx);

                // if:
                // check 3D only
                // check with the last component element computed for the current 3D element
                // intersection has not been computed already
                if (ele_3D->dim() == 3 &&
                    (last_slave_for_3D_elements[bulk_ele_idx] != component_ele_idx &&
                     !intersection_exists(component_ele_idx,bulk_ele_idx) )
                ) {
                    // check that tetrahedron element is numbered correctly and is not degenerated
                    ASSERT_DBG(ele_3D->tetrahedron_jacobian() > 0).add_value(ele_3D->index(),"element index").error(
                           "Tetrahedron element (%d) has wrong numbering or is degenerated (negative Jacobian).");
                    
                        // - find first intersection
                        // - if found, prolongate and possibly fill both prolongation queues
                        // do-while loop:
                        // - empty prolongation queues:
                        //      - empty bulk queue:
                        //          - get a candidate from queue and compute CI
                        //          - prolongate and possibly push new candidates into queues
                        //          - repeat until bulk queue is empty
                        //          - the component element is still the same whole time in here
                        //
                        //      - the component element might get fully covered by bulk elements
                        //        and only then it can be closed
                        //
                        //      - pop next candidate from component queue:
                        //          - the component element is now changed
                        //          - compute CI
                        //          - prolongate and possibly push new candidates into queues
                        //
                        // - repeat until both queues are empty
                    
                    update_simplex(elm, component_simplex); // update component simplex
                    update_simplex(ele_3D, tetrahedron); // update tetrahedron
                    std::vector<unsigned int> prolongation_table;
                    bool found = compute_initial_CI(component_ele_idx, bulk_ele_idx,
                                                    component_counter_, prolongation_table);

                    // keep the index of the current component element that is being investigated
                    unsigned int current_component_element_idx = component_ele_idx;
                    
                    if(found){
                        component_counter_++;
                        DebugOut().fmt("comp: {}\n", component_counter_);
                        
                        DebugOut().fmt("start component with elements {} {}\n",component_ele_idx, bulk_ele_idx);
                        
                        prolongation_decide(elm, ele_3D, intersection_list_[component_ele_idx].back(), prolongation_table);
                        
                        START_TIMER("Prolongation algorithm");
                        do{
                            // flag is set false if the component element is not fully covered with tetrahedrons
                            bool element_covered = true;
                            
                            while(!bulk_queue_.empty()){
                                Prolongation pr = bulk_queue_.front();
                                DebugOut().fmt("Bulk queue: ele_idx {}.\n",pr.elm_3D_idx);
                                
                                if( pr.elm_3D_idx == undefined_elm_idx_)
                                {
                                    DebugOut().fmt("Open intersection component element: {}\n",current_component_element_idx);
                                    element_covered = false;
                                }
                                else prolongate(pr);
                                
                                bulk_queue_.pop();
                            }
                            
                            if(! closed_elements[current_component_element_idx])
                                closed_elements[current_component_element_idx] = element_covered;
                            
                            
                            if(!component_queue_.empty()){
                                Prolongation pr = component_queue_.front();

                                // note the component element index
                                current_component_element_idx = pr.component_elm_idx;
                                DebugOut().fmt("Component queue: ele_idx {}.\n",current_component_element_idx);
                                
                                prolongate(pr);
                                component_queue_.pop();
                            }
                        }
                        while( !(component_queue_.empty() && bulk_queue_.empty()) );
                        END_TIMER("Prolongation algorithm");
                        
                        // if component element is closed, do not check other bounding boxes
                        if(closed_elements[component_ele_idx])
                            break;
                    }

                }
            }
            END_TIMER("Bounding box element iteration");
        }
    }

    END_TIMER("Element iteration");
    
    // DBG write which elements are closed
//     FOR_ELEMENTS(mesh, ele) {
//         DebugOut().fmt("Element[{}] closed: {}\n",ele.index(),(closed_elements[ele.index()] ? 1 : 0));
//     }
}
  
template<unsigned int dim>
void InspectElementsAlgorithm<dim>::compute_intersections_BIHtree(const BIHTree& bih)
{
    DebugOut() << "#########   ALGORITHM: compute_intersections_BIHtree   #########\n";
    
    init();
    
    START_TIMER("Element iteration");
    
    FOR_ELEMENTS(mesh, elm) {
        unsigned int component_ele_idx = elm->index();
        
        if (elm->dim() == dim &&                                    // is component element
            bih.ele_bounding_box(component_ele_idx).intersect(bih.tree_box()))   // its bounding box intersects 3D mesh bounding box
        {   
            update_simplex(elm, component_simplex); // update component simplex
            std::vector<unsigned int> searchedElements;
            
            START_TIMER("BIHtree find");
            bih.find_bounding_box(bih.ele_bounding_box(component_ele_idx), searchedElements);
            END_TIMER("BIHtree find");
            
            START_TIMER("Bounding box element iteration");
            
            // Go through all element which bounding box intersects the component element bounding box
            for (std::vector<unsigned int>::iterator it = searchedElements.begin(); it!=searchedElements.end(); it++)
            {
                unsigned int bulk_ele_idx = *it;
                ElementFullIter ele_3D = mesh->element(bulk_ele_idx);
                
                if (ele_3D->dim() == 3
                ) {
                    // check that tetrahedron element is numbered correctly and is not degenerated
                    ASSERT_DBG(ele_3D->tetrahedron_jacobian() > 0).add_value(ele_3D->index(),"element index").error(
                           "Tetrahedron element (%d) has wrong numbering or is degenerated (negative Jacobian).");
                    
                    update_simplex(ele_3D, tetrahedron); // update tetrahedron
                    std::vector<unsigned int> prolongation_table;
                    
                    IntersectionAux<dim,3> is(component_ele_idx, bulk_ele_idx, 0);
                    START_TIMER("Compute intersection");
                    ComputeIntersection<Simplex<dim>, Simplex<3>> CI(component_simplex, tetrahedron);
                    CI.init();
                    CI.compute(is, prolongation_table);
                    END_TIMER("Compute intersection");
                    
                    if(is.points().size() > 0) {
                        
                        intersection_list_[component_ele_idx].push_back(is);
                        n_intersections_++;
                        // if component element is closed, do not check other bounding boxes
                        closed_elements[component_ele_idx] = true;
                    }
                }
            }
            END_TIMER("Bounding box element iteration");
        }
    }

    END_TIMER("Element iteration");
}

template<unsigned int dim>
void InspectElementsAlgorithm<dim>::compute_intersections_BB()
{
    DebugOut() << "#########   ALGORITHM: compute_intersections_BB   #########\n";
    init();
    compute_bounding_boxes();
    
    START_TIMER("Element iteration");
    
    
    FOR_ELEMENTS(mesh, elm) {
        unsigned int component_ele_idx = elm->index();
        
        if (elm->dim() == dim &&                                // is component element
            !closed_elements[component_ele_idx] &&                    // is not closed yet
            elements_bb[component_ele_idx].intersect(mesh_3D_bb))    // its bounding box intersects 3D mesh bounding box
        {    
            
            START_TIMER("Bounding box element iteration");
            
            // Go through all element which bounding box intersects the component element bounding box
            FOR_ELEMENTS(mesh, ele_3D)
            {
                unsigned int bulk_ele_idx = ele_3D.index();

                // if:
                // check 3D only
                // check with the last component element computed for the current 3D element
                // check that the bounding boxes intersect
                // intersection has not been computed already
                if (ele_3D->dim() == 3 &&
                    (last_slave_for_3D_elements[bulk_ele_idx] != component_ele_idx &&
                     elements_bb[component_ele_idx].intersect(elements_bb[bulk_ele_idx]) &&
                     !intersection_exists(component_ele_idx,bulk_ele_idx) )
                ){
                    // check that tetrahedron element is numbered correctly and is not degenerated
                    ASSERT_DBG(ele_3D->tetrahedron_jacobian() > 0).add_value(ele_3D->index(),"element index").error(
                           "Tetrahedron element (%d) has wrong numbering or is degenerated (negative Jacobian).");
                    
                        // - find first intersection
                        // - if found, prolongate and possibly fill both prolongation queues
                        // do-while loop:
                        // - empty prolongation queues:
                        //      - empty bulk queue:
                        //          - get a candidate from queue and compute CI
                        //          - prolongate and possibly push new candidates into queues
                        //          - repeat until bulk queue is empty
                        //          - the component element is still the same whole time in here
                        //
                        //      - the component element might get fully covered by bulk elements
                        //        and only then it can be closed
                        //
                        //      - pop next candidate from component queue:
                        //          - the component element is now changed
                        //          - compute CI
                        //          - prolongate and possibly push new candidates into queues
                        //
                        // - repeat until both queues are empty
                    
                    update_simplex(elm, component_simplex); // update component simplex
                    update_simplex(ele_3D, tetrahedron); // update tetrahedron
                    std::vector<unsigned int> prolongation_table;
                    bool found = compute_initial_CI(component_ele_idx, bulk_ele_idx,
                                                    component_counter_, prolongation_table);

                    // keep the index of the current component element that is being investigated
                    unsigned int current_component_element_idx = component_ele_idx;
                    
                    if(found){
                        component_counter_++;
                        DebugOut().fmt("comp: {}\n", component_counter_);
                        
                        DebugOut().fmt("start component with elements {} {}\n",component_ele_idx, bulk_ele_idx);
                        
                        prolongation_decide(elm, ele_3D, intersection_list_[component_ele_idx].back(), prolongation_table);
                        
                        START_TIMER("Prolongation algorithm");
                        do{
                            // flag is set false if the component element is not fully covered with tetrahedrons
                            bool element_covered = true;
                            
                            while(!bulk_queue_.empty()){
                                Prolongation pr = bulk_queue_.front();
                                DebugOut().fmt("Bulk queue: ele_idx {}.\n",pr.elm_3D_idx);
                                
                                if( pr.elm_3D_idx == undefined_elm_idx_)
                                {
                                    DebugOut().fmt("Open intersection component element: {}\n",current_component_element_idx);
                                    element_covered = false;
                                }
                                else prolongate(pr);
                                
                                bulk_queue_.pop();
                            }
                            
                            if(! closed_elements[current_component_element_idx])
                                closed_elements[current_component_element_idx] = element_covered;
                            
                            
                            if(!component_queue_.empty()){
                                Prolongation pr = component_queue_.front();

                                // note the component element index
                                current_component_element_idx = pr.component_elm_idx;
                                DebugOut().fmt("Component queue: ele_idx {}.\n",current_component_element_idx);
                                
                                prolongate(pr);
                                component_queue_.pop();
                            }
                        }
                        while( !(component_queue_.empty() && bulk_queue_.empty()) );
                        END_TIMER("Prolongation algorithm");
                        
                        // if component element is closed, do not check other bounding boxes
                        if(closed_elements[component_ele_idx])
                            break;
                    }

                }
            }
            END_TIMER("Bounding box element iteration");
        }
    }

    END_TIMER("Element iteration");
    
    // DBG write which elements are closed
//     FOR_ELEMENTS(mesh, ele) {
//         DebugOut().fmt("Element[{}] closed: {}\n",ele.index(),(closed_elements[ele.index()] ? 1 : 0));
//     }
}



template<unsigned int dim>
template<unsigned int ele_dim>
std::vector< unsigned int > InspectElementsAlgorithm<dim>::get_element_edges(const ElementFullIter& ele,
                                                                             unsigned int ip_dim,
                                                                             unsigned int ip_obj_idx,
                                                                             const bool &include_current_ele
                                                                             )
{
    std::vector<Edge*> edges;
    edges.reserve(ele_dim - ip_dim);  // reserve number of possible edges

    DebugOut() << "dim " << ele_dim << ": ";
    switch (ip_dim)
    {
        // IP is at a node of tetrahedron; possible edges are from all connected sides (3)
        case 0: if(ele_dim == 1) {
                    edges.push_back(&(mesh->edges[ele->edge_idx_[ip_obj_idx]]));
                    break;
                }
                
                for(unsigned int j=0; j < RefElement<ele_dim>::n_sides_per_node; j++){
                    unsigned int local_edge = RefElement<ele_dim>::interact(Interaction<ele_dim-1,0>(ip_obj_idx))[j];
                    edges.push_back(&(mesh->edges[ele->edge_idx_[local_edge]]));
                }
                DebugOut() << "prolong (node)\n";
                break;
        
        // IP is on a line of tetrahedron; possible edges are from all connected sides (2)
        case 1: if(ele_dim == 2) {
                    edges.push_back(&(mesh->edges[ele->edge_idx_[ip_obj_idx]]));
                    break;
                }
            
                ASSERT_DBG(ele_dim == 3);
                for(unsigned int j=0; j < RefElement<ele_dim>::n_sides_per_line; j++){
                    unsigned int local_edge = RefElement<ele_dim>::interact(Interaction<2,1>(ip_obj_idx))[j];
                    edges.push_back(&(mesh->edges[ele->edge_idx_[local_edge]]));
                }
                DebugOut() << "prolong (edge)\n";
                break;
                
        // IP is on a side of tetrahedron; only possible edge is from the given side (1)
        case 2: ASSERT_DBG(ele_dim == 3);
                edges.push_back(&(mesh->edges[ele->edge_idx_[ip_obj_idx]]));
                DebugOut() << "prolong (side)\n";
                break;
        default: ASSERT_DBG(0);
    }
    
    // get indices of neighboring bulk elements
    std::vector<unsigned int> elements_idx;
    elements_idx.reserve(2*edges.size());    // twice the number of edges
    for(Edge* edg : edges)
    for(int j=0; j < edg->n_sides;j++) {
        if (edg->side(j)->element() != ele)
            elements_idx.push_back(edg->side(j)->element()->index());
    }
    
    // possibly include the current bulk element
    if(include_current_ele)
        elements_idx.push_back(ele->index());
    
    return elements_idx;
}

template<unsigned int dim>
unsigned int InspectElementsAlgorithm<dim>::create_prolongations_over_bulk_element_edges(const std::vector< unsigned int >& bulk_neighbors,
                                                                                         const unsigned int &component_ele_idx)
{
    unsigned int n_prolongations = 0;
    for(unsigned int bulk_neighbor_idx : bulk_neighbors)
    {
        if(last_slave_for_3D_elements[bulk_neighbor_idx] == undefined_elm_idx_ ||
            (last_slave_for_3D_elements[bulk_neighbor_idx] != component_ele_idx && !intersection_exists(component_ele_idx,bulk_neighbor_idx)))
        {
            last_slave_for_3D_elements[bulk_neighbor_idx] = component_ele_idx;
        
            DebugOut().fmt("3d prolong {} in {}\n",component_ele_idx,bulk_neighbor_idx);
            
            // Vytvoření průniku bez potřeby počítání
            IntersectionAux<dim,3> il_other(component_ele_idx, bulk_neighbor_idx);
            intersection_list_[component_ele_idx].push_back(il_other);
            
            Prolongation pr = {component_ele_idx, bulk_neighbor_idx, (unsigned int)intersection_list_[component_ele_idx].size() - 1};
            bulk_queue_.push(pr);
            n_prolongations++;
        }
    }
    return n_prolongations;
}

// template<unsigned int dim>
// unsigned int InspectElementsAlgorithm<dim>::create_prolongations_over_bulk_element_edges(const unsigned int& bulk_ele_idx,
//                                                                                          const unsigned int& component_ele_idx)
// {
//     if(last_slave_for_3D_elements[bulk_ele_idx] == undefined_elm_idx_ ||
//         (last_slave_for_3D_elements[bulk_ele_idx] != component_ele_idx && !intersection_exists(component_ele_idx,bulk_ele_idx)))
//     {
//         last_slave_for_3D_elements[bulk_ele_idx] = component_ele_idx;
//     
//         DebugOut().fmt("bulk prolongation: c {} in b {}\n",component_ele_idx,bulk_ele_idx);
//         
//         // Vytvoření průniku bez potřeby počítání
//         IntersectionAux<dim,3> il_other(component_ele_idx, bulk_ele_idx);
//         intersection_list_[component_ele_idx].push_back(il_other);
//         
//         Prolongation pr = {component_ele_idx, bulk_ele_idx, (unsigned int)intersection_list_[component_ele_idx].size() - 1};
//         bulk_queue_.push(pr);
//         
//         return 1;
//     }
//     return 0;
// }


template<>
void InspectElementsAlgorithm<1>::prolongation_decide(const ElementFullIter& comp_ele,
                                                      const ElementFullIter& bulk_ele,
                                                      const IntersectionAux<1,3> &is,
                                                      const std::vector<unsigned int> &prolongation_table)
{
    DebugOut() << "DECIDE\n";
    // number of IPs that are at vertices of component element (counter used for closing element)
    unsigned int n_ip_vertices = 0;
    
    for(const IntersectionPointAux<1,3> &IP : is.points()) {
        if(IP.dim_A() == 0) { // if IP is the end of the 1D element
            n_ip_vertices++;
            DebugOut() << "1D end\n";
            
            // there are two possibilities:
            // 1] IP is inside the bulk element
            //    => prolongate 1D element as long as it creates prolongation point on the side of tetrahedron
            // 2] IP lies on the boundary of tetrahedron (vertex, edge or side)
            //    => search all connected edges of tetrahedron for neghboring sides
            //       and create candidates: neighboring 1D element + neighboring tetrahedron
                
            if(IP.dim_B() == 3)
            {
                // iterate over sides of 1D element
                SideIter elm_side = comp_ele->side(IP.idx_A());  // side of 1D element is vertex
                Edge *edg = elm_side->edge();
                for(int j=0; j < edg->n_sides;j++) {
                    
                    SideIter other_side=edg->side(j);
                    if (other_side != elm_side) {   // we do not want to look at the current 1D element

                        unsigned int component_neighbor_idx = other_side->element()->index();
                        if(!intersection_exists(component_neighbor_idx,bulk_ele->index())){
                            DebugOut().fmt("1d prolong {} in {}\n", component_neighbor_idx, bulk_ele->index());
                            
                                // Vytvoření průniku bez potřeby počítání
                                IntersectionAux<1,3> il_other(component_neighbor_idx, bulk_ele->index());
                                intersection_list_[component_neighbor_idx].push_back(il_other);
                                
                                Prolongation pr = {component_neighbor_idx, bulk_ele->index(), 
                                                (unsigned int)intersection_list_[component_neighbor_idx].size() - 1};
                                component_queue_.push(pr);
                        }
                    }
                }
            }
            else
            {
                // search for indices of neighboring bulk elements (including the current one)
                std::vector<unsigned int> bulk_neighbors = get_element_edges<3>(bulk_ele,IP.dim_B(), IP.idx_B(), true);
                
                unsigned int n_prolongations = 0;
                // iterate over sides of 1D element
                SideIter elm_side = comp_ele->side(IP.idx_A());  // side of 1D element is vertex
                Edge *edg = elm_side->edge();
                for(int j=0; j < edg->n_sides;j++) {
                    
                    SideIter side=edg->side(j);
                    unsigned int component_neighbor_idx = side->element()->index();
                    // we want to look also at the current 1D and 3D element
                    
                    n_prolongations += create_prolongations_over_bulk_element_edges(bulk_neighbors,component_neighbor_idx);
                }
                
                // if there are no sides of any edge that we can continue to prolongate over,
                // it means we are at the boundary and cannot prolongate further
                if(n_prolongations == 0)
                {
                    Prolongation pr = {comp_ele->index(), undefined_elm_idx_, undefined_elm_idx_};
                    bulk_queue_.push(pr);
                }
            }

        }else{
            std::vector<unsigned int> bulk_neighbors = get_element_edges<3>(bulk_ele,IP.dim_B(), IP.idx_B(),false);
            
            // if edge has only one side, it means it is on the boundary and we cannot prolongate
            if(bulk_neighbors.empty())
            {
                Prolongation pr = {comp_ele->index(), undefined_elm_idx_, undefined_elm_idx_};
                bulk_queue_.push(pr);
                continue;
            }

            unsigned int n_prolongations = create_prolongations_over_bulk_element_edges(bulk_neighbors,comp_ele->index());
            
            DebugOut().fmt("cover: {} {}\n", is.size(), n_prolongations);
            // if there are no sides of any edge that we can continue to prolongate over,
            // it means we are at the boundary and cannot prolongate further
            if(bulk_neighbors.size() != 1 && n_prolongations == 0)
            {
                Prolongation pr = {comp_ele->index(), undefined_elm_idx_, undefined_elm_idx_};
                bulk_queue_.push(pr);
            }
        }  
    }
    
    // close component element if it has all vertices inside bulk element
    if(n_ip_vertices == is.size()) closed_elements[comp_ele->index()] = true;
}


template<>
void InspectElementsAlgorithm<2>::assert_same_intersection(unsigned int comp_ele_idx, unsigned int bulk_ele_idx)
{
    for(unsigned int i=0; i < intersection_list_[comp_ele_idx].size(); i++)
    {
        if(intersection_list_[comp_ele_idx][i].bulk_ele_idx() == bulk_ele_idx)
        {
            DebugOut().fmt("intersection comp-bulk: {} {}\n", comp_ele_idx, bulk_ele_idx);
            ASSERT_DBG(0).add_value(bulk_ele_idx,"bulk_ele_idx").error("Want to add the same intersection!");
        }
    }
}

template<>
void InspectElementsAlgorithm<2>::prolongation_decide(const ElementFullIter& comp_ele,
                                                      const ElementFullIter& bulk_ele,
                                                      const IntersectionAux<2,3> &is,
                                                      const std::vector<unsigned int> &prolongation_table)
{
    DebugOut() << "DECIDE\n";
    // number of IPs that are at vertices of component element (counter used for closing element)
    unsigned int n_ip_vertices = 0;
    
    for(const IntersectionPointAux<2,3> &IP : is.points()) {
        if(IP.dim_A() < 2) { // if IP on one of the sides of triangle
            if(IP.dim_A() == 0) n_ip_vertices++;
            DebugOut() << "on 2D side, dim = " << IP.dim_A() << "\n";
            
            // search for indices of neighboring component elements (including the current one)
            std::vector<unsigned int> comp_neighbors = get_element_edges<2>(comp_ele,IP.dim_A(), IP.idx_A(),false);
                
            // there are two possibilities:
            // 1] IP is inside the bulk element
            //    => prolongate 1D element as long as it creates prolongation point on the side of tetrahedron
            // 2] IP lies on the boundary of tetrahedron (vertex, edge or side)
            //    => search all connected edges of tetrahedron for neghboring sides
            //       and create candidates: neighboring 1D element + neighboring tetrahedron
                
            if(IP.dim_B() == 3)
            {
                for(unsigned int j=0; j < comp_neighbors.size(); j++) {

                    unsigned int comp_neighbor_idx = comp_neighbors[j],
                                 bulk_current = bulk_ele->index();
                    if(!intersection_exists(comp_neighbor_idx,bulk_current)){
                        DebugOut().fmt("2d prolong {} in {}\n", comp_neighbor_idx, bulk_current);
                        
                            // Vytvoření průniku bez potřeby počítání
                            IntersectionAux<2,3> il_other(comp_neighbor_idx, bulk_current);
                            intersection_list_[comp_neighbor_idx].push_back(il_other);
                            
                            Prolongation pr = {comp_neighbor_idx, bulk_current, 
                                            (unsigned int)intersection_list_[comp_neighbor_idx].size() - 1};
                            component_queue_.push(pr);
                            last_slave_for_3D_elements[bulk_current] = comp_neighbor_idx;
                    }
                }
            }
            else
            {
                // search for indices of neighboring bulk elements (including the current one)
                std::vector<unsigned int> bulk_neighbors = get_element_edges<3>(bulk_ele,IP.dim_B(),IP.idx_B(), true);
                
                unsigned int n_prolongations = 0;
                // prolong also current element
                n_prolongations += create_prolongations_over_bulk_element_edges(bulk_neighbors,comp_ele->index());
                
                for(unsigned int j=0; j < comp_neighbors.size(); j++) {
                    n_prolongations += create_prolongations_over_bulk_element_edges(bulk_neighbors,comp_neighbors[j]);
                }
                
                // if there are no sides of any edge that we can continue to prolongate over,
                // it means we are at the boundary and cannot prolongate further
                if(n_prolongations == 0)
                {
                    Prolongation pr = {comp_ele->index(), undefined_elm_idx_, undefined_elm_idx_};
                    bulk_queue_.push(pr);
                }
            }

        }else{
            std::vector<unsigned int> bulk_neighbors = get_element_edges<3>(bulk_ele,IP.dim_B(), IP.idx_B(),false);
            
            // if edge has only one side, it means it is on the boundary and we cannot prolongate
            if(bulk_neighbors.empty())
            {
                Prolongation pr = {comp_ele->index(), undefined_elm_idx_, undefined_elm_idx_};
                bulk_queue_.push(pr);
                continue;
            }

            unsigned int n_prolongations = create_prolongations_over_bulk_element_edges(bulk_neighbors,comp_ele->index());
            
            DebugOut().fmt("cover: {} {}\n", is.size(), n_prolongations);
            // if there are no sides of any edge that we can continue to prolongate over,
            // it means we are at the boundary and cannot prolongate further
            if(bulk_neighbors.size() != 1 && n_prolongations == 0)
            {
                Prolongation pr = {comp_ele->index(), undefined_elm_idx_, undefined_elm_idx_};
                bulk_queue_.push(pr);
            }
        }  
    }
    
    // close component element if it has all vertices inside bulk element
    if(n_ip_vertices == is.size()) closed_elements[comp_ele->index()] = true;
}

// template<>
// void InspectElementsAlgorithm<2>::prolongation_decide(const ElementFullIter& comp_ele,
//                                                       const ElementFullIter& bulk_ele,
//                                                       const IntersectionAux<2,3> &is,
//                                                       const std::vector<unsigned int> &prolongation_table
//                                                      )
// {
//     DebugOut() << "DECIDE\n";
//     // in case there are more than 2 IPs, the prolongation table is filled
//     for(unsigned int i = 0; i < prolongation_table.size();i++){
// 
//         unsigned int side;
//         bool is_triangle_side = true;
// 
//         if(prolongation_table[i] >= 4){
//             side = prolongation_table[i] - 4;
//         }else{
//             side = prolongation_table[i];
//             is_triangle_side = false;
//         }
// 
//         DebugOut().fmt("prolongation table: {} {}\n", side, is_triangle_side);
// 
//         if(is_triangle_side){
//             // prolongation through the triangle side
// 
//             SideIter elm_side = comp_ele->side(side);
//             Edge *edg = elm_side->edge();
// 
//             for(int j=0; j < edg->n_sides;j++) {
//                 SideIter other_side=edg->side(j);
//                 if (other_side != elm_side) {
//                     unsigned int sousedni_element = other_side->element()->index(); // 2D element
// 
//                     DebugOut().fmt("2d sousedni_element {}\n", sousedni_element);
// 
//                         if(!intersection_exists(sousedni_element,bulk_ele->index())){
// 
//                             unsigned int component_idx = is.component_idx();
//                             if(edg->n_sides > 2)
//                             {
//                                 component_counter_++;
//                                 DebugOut().fmt("COMPONENT COUNTER: {}\n", component_counter_);
//                                 component_idx = component_counter_;
//                             }
//                             
//                             DebugOut() << "2d prolong\n";
//                             assert_same_intersection(sousedni_element, bulk_ele->index());
//                             // Vytvoření průniku bez potřeby počítání
//                             IntersectionAux<2,3> il_other(sousedni_element, bulk_ele->index(), component_idx);
//                             intersection_list_[sousedni_element].push_back(il_other);
// 
//                             Prolongation pr = {sousedni_element, bulk_ele->index(), (unsigned int)intersection_list_[sousedni_element].size() - 1};
//                             component_queue_.push(pr);
//                         }
//                 }
//             }
// 
//         }else{
//             // prolongation through the tetrahedron side
//             
//             SideIter elm_side = bulk_ele->side(side);
//             Edge *edg = elm_side->edge();
// 
//             for(int j=0; j < edg->n_sides;j++) {
//                 SideIter other_side=edg->side(j);
//                 if (other_side != elm_side) {
//                     //
// 
//                     unsigned int sousedni_element = other_side->element()->index();
// 
//                     DebugOut().fmt("3d sousedni_element {}\n", sousedni_element);
// 
//                     // TODO:
//                     // - rename it, describe it and test that it is really useful !!
//                     // how it probably works: 
//                     // - if index undefined then no intersection has been computed for the 3D element
//                     // - if not undefined check the index of actual 2D element 
//                     //   (most probable case, that we are looking along 2D element back to 3D element, which we have just computed)
//                     // - if it is another 2D element, then go through all found intersections of the 3D element and test it..
//                     
//                     if(last_slave_for_3D_elements[sousedni_element] == undefined_elm_idx_ || 
//                         (last_slave_for_3D_elements[sousedni_element] != comp_ele->index() && !intersection_exists(comp_ele->index(),sousedni_element))){
//                         
//                         last_slave_for_3D_elements[sousedni_element] = comp_ele->index();
// 
//                         DebugOut() << "3d prolong\n";
//                         DebugOut().fmt("last_slave_for_3D_elements: {}, exists: {}\n", last_slave_for_3D_elements[sousedni_element], intersection_exists(comp_ele->index(),sousedni_element));
//                         assert_same_intersection(comp_ele->index(), sousedni_element);
//                         // Vytvoření průniku bez potřeby počítání
//                         IntersectionAux<2,3> il_other(comp_ele->index(), sousedni_element, is.component_idx());
//                         intersection_list_[comp_ele->index()].push_back(il_other);
// 
//                         Prolongation pr = {comp_ele->index(), sousedni_element, (unsigned int)intersection_list_[comp_ele->index()].size() - 1};
//                         bulk_queue_.push(pr);
// 
//                     }
//                 }
//             }
//         }
//     }
//     // otherwise, we need to be able to prolongate over 1 or 2 IPs; if we want to find components properly
//     // TODO: ...
// }

template<unsigned int dim>
void InspectElementsAlgorithm<dim>::prolongate(const InspectElementsAlgorithm< dim >::Prolongation& pr)
{
    ElementFullIter elm = mesh->element(pr.component_elm_idx);
    ElementFullIter ele_3D = mesh->element(pr.elm_3D_idx);
    
    DebugOut().fmt("Prolongate {}D: {} in {}.\n", dim, pr.component_elm_idx, pr.elm_3D_idx);

    //TODO: optimization: this might be called before and not every time 
    //(component element is not changing when emptying bulk queue)
    update_simplex(elm, component_simplex);
    update_simplex(ele_3D, tetrahedron);

    IntersectionAux<dim,3> &is = intersection_list_[pr.component_elm_idx][pr.dictionary_idx];
    std::vector<unsigned int> prolongation_table;
    
    START_TIMER("Compute intersection");
    ComputeIntersection<Simplex<dim>, Simplex<3>> CI(component_simplex, tetrahedron);
    CI.init();
    CI.compute(is, prolongation_table);
    END_TIMER("Compute intersection");
    
    last_slave_for_3D_elements[pr.elm_3D_idx] = pr.component_elm_idx;
    
    if(is.size() > 0){
//         for(unsigned int j=0; j < is.size(); j++) 
//             DebugOut() << is[j];
        
        prolongation_decide(elm, ele_3D,is,prolongation_table);
        n_intersections_++;
    }
}




InspectElementsAlgorithm22::InspectElementsAlgorithm22(Mesh* input_mesh)
: mesh(input_mesh)
{}


void InspectElementsAlgorithm22::compute_intersections(const std::vector< std::vector<ILpair>>& intersection_map_)
{
    DebugOut() << "Intersections 2d-2d\n";
    
    FOR_ELEMENTS(mesh, ele) {
    if (ele->dim() == 3)
    {
        unsigned int ele_idx = ele->index();
        // if there are not at least 2 2D elements intersecting 3D element; continue
        if(intersection_map_[ele_idx].size() < 2) continue;
        
        const std::vector<ILpair> &local_map = intersection_map_[ele_idx];
        
        //DebugOut() << "more than 2 intersections in tetrahedron found\n";
        for(unsigned int i=0; i < local_map.size(); i++)
        {
            //TODO: 1] compute all plucker coords at once
            //TODO: 2] pass plucker coords from 2d-3d
            
            ElementFullIter eleA = mesh->element(local_map[i].first);
            if(eleA->dim() !=2 ) continue;  //skip other dimension intersection
            unsigned int componentA_idx = local_map[i].second->component_idx();
            
            IntersectionLocalBase * ilb = local_map[i].second;
            DebugOut().fmt("2d-2d ILB: {} {} {}\n", ilb->bulk_ele_idx(), ilb->component_ele_idx(), ilb->component_idx());
            
            for(unsigned int j=i+1; j < local_map.size(); j++)
            {
                ElementFullIter eleB = mesh->element(local_map[j].first);
                if(eleB->dim() !=2 ) continue;  //skip other dimension intersection
                
                // component check not working, until prolongation will be done also over vertices..
                unsigned int componentB_idx = local_map[j].second->component_idx();
                if(componentA_idx == componentB_idx) continue;  //skip elements of the same component
                // this also skips the compatible connections (it is still a single component in this case)
                
//                 //does not solve 'vertex neighbors' (common only one node)
//                 bool is_not_neighbor = true;
//                 for(unsigned int k=0; k < eleA->n_sides(); k++)
//                 {
//                     Edge * edge = eleA->side(k)->edge();
//                     for(unsigned int s=0; s < edge->n_sides; s++)
//                     {
//                         if(eleA->side(k) != edge->side(s))
//                             if(edge->side(s)->element()->index() == eleB.index()) is_not_neighbor = false;
//                     }
//                 }
//                 if(is_not_neighbor) continue;
                
                DebugOut().fmt("compute intersection 2d-2d: e_{} e_{} c_{} c_{}\n",eleA.index(), eleB.index(), componentA_idx, componentB_idx);
                compute_single_intersection(eleA,
                                            eleB);
            }
        }
    }
    }
}

void InspectElementsAlgorithm22::compute_single_intersection(const ElementFullIter& eleA,
                                                             const ElementFullIter& eleB)
{
    ASSERT_DBG(eleA->dim() == 2);
    ASSERT_DBG(eleB->dim() == 2);
    ASSERT_DBG(eleA->index() != eleB->index());
    
    update_simplex(eleA, triaA_);
    update_simplex(eleB, triaB_);
    
    IntersectionAux<2,2> is(eleA->index(), eleB->index(), 0);
    std::vector<unsigned int> prolongation_table;
    
    ComputeIntersection< Simplex<2>, Simplex<2>> CI(triaA_, triaB_);
    CI.init();
    unsigned int n_local_intersection = CI.compute(is, prolongation_table);
    
    if(n_local_intersection > 0)
        intersectionaux_storage22_.push_back(is);
}





InspectElementsAlgorithm12::InspectElementsAlgorithm12(Mesh* input_mesh)
: mesh(input_mesh)
{}


void InspectElementsAlgorithm12::compute_intersections(std::vector< std::vector<ILpair>>& intersection_map,
                                                       std::vector<IntersectionLocal<1,2>> &storage)
{
    DebugOut() << "Intersections 1d-2d\n";
    
    FOR_ELEMENTS(mesh, ele) {
    if (ele->dim() == 3)
    {
        unsigned int ele_idx = ele->index();
        // if there are not at least 2 elements intersecting 3D element; continue
        if(intersection_map[ele_idx].size() < 2) continue;
        
        const std::vector<ILpair> &local_map = intersection_map[ele_idx];
        
        //DebugOut() << "more than 2 intersections in tetrahedron found\n";
        for(unsigned int i=0; i < local_map.size(); i++)
        {
            //TODO: 1] compute all plucker coords at once
            //TODO: 2] pass plucker coords from 1d-3d
            
            unsigned int eleA_idx = local_map[i].first;
            ElementFullIter eleA = mesh->element(eleA_idx);
            
            if(eleA->dim() !=1 ) continue;  //skip other dimension intersection
            
            for(unsigned int j=0; j < local_map.size(); j++)
            {
                unsigned int eleB_idx = local_map[j].first;
                ElementFullIter eleB = mesh->element(local_map[j].first);
                if(eleB->dim() !=2 ) continue;  //skip other dimension intersection
                
                //skip candidates already computed
                bool skip = false;
                for(unsigned int i=0; i<intersection_map[eleA_idx].size(); i++)
                {
                    if(intersection_map[eleA_idx][i].first == eleB_idx) {
                        skip = true;
                        break;
                    }
                }
                
                if(skip) continue;
                
                DebugOut().fmt("compute intersection 1d-2d: {} {}\n",eleA.index(), eleB.index());
//                 compute_single_intersection(eleA,
//                                             eleB);
                
                update_simplex(eleA, abscissa_);
                update_simplex(eleB, triangle_);
                
                IntersectionAux<1,2> is(eleA_idx, eleB_idx, 0);
                
                ComputeIntersection< Simplex<1>, Simplex<2>> CI(abscissa_, triangle_);
                unsigned int n_local_intersection = CI.compute_final(is.points());
    
                if(n_local_intersection > 0)
                {
                    storage.push_back(IntersectionLocal<1,2>(is));
                    intersection_map[eleA_idx].push_back(std::make_pair(
                                                         eleB_idx,
                                                         &(storage.back())
                                                         ));
                    intersection_map[eleB_idx].push_back(std::make_pair(
                                                         eleA_idx,
                                                         &(storage.back())
                                                         ));
                    
                    DebugOut().fmt("1D-2D intersection [{} - {}]:\n",is.component_ele_idx(), is.bulk_ele_idx());
                }
            }
        }
    }
    }
    
    // just dbg output
    for(IntersectionLocal<1,2> &is : storage)
    {
        DebugOut().fmt("1D-2D intersection [{} - {}]:\n",is.component_ele_idx(), is.bulk_ele_idx());
        for(const IntersectionPoint<1,2>& ip : is.points()) {
            //DebugOut() << ip;
            auto p = ip.coords(mesh->element(is.component_ele_idx()));
            DebugOut() << "[" << p[0] << " " << p[1] << " " << p[2] << "]\n";
        }
    }
}

// void InspectElementsAlgorithm12::compute_single_intersection(const ElementFullIter& eleA,
//                                                              const ElementFullIter& eleB)
// {
//     ASSERT_DBG(eleA->dim() == 1);
//     ASSERT_DBG(eleB->dim() == 2);
//     
//     update_simplex(eleA, abscissa_);
//     update_simplex(eleB, triangle_);
//     
//     IntersectionAux<1,2> is(eleA->index(), eleB->index(), 0);
// //     std::vector<unsigned int> prolongation_table;
//     
//     ComputeIntersection< Simplex<1>, Simplex<2>> CI(abscissa_, triangle_);
//     unsigned int n_local_intersection = CI.compute_final(is.points());
//     
//     if(n_local_intersection > 0)
//     {
//         DebugOut().fmt("found: {}\n",n_local_intersection);
//         intersectionaux_storage12_.push_back(is);
//     }
// }



void InspectElementsAlgorithm12::compute_intersections_2(const BIHTree& bih)
{
    DebugOut() << "Intersections 1d-2d (2-bihtree)\n";
    
    START_TIMER("Element iteration");
    
    FOR_ELEMENTS(mesh, elm) {
        unsigned int component_ele_idx = elm->index();
        
        if (elm->dim() == 1)                                    // is component element
            //&& elements_bb[component_ele_idx].intersect(mesh_3D_bb))   // its bounding box intersects 3D mesh bounding box
        {   
            update_simplex(elm, abscissa_); // update component simplex
            std::vector<unsigned int> searchedElements;
            
            START_TIMER("BIHtree find");
            bih.find_bounding_box(bih.ele_bounding_box(component_ele_idx), searchedElements);
            END_TIMER("BIHtree find");
            
            START_TIMER("Bounding box element iteration");
            
            // Go through all element which bounding box intersects the component element bounding box
            for (std::vector<unsigned int>::iterator it = searchedElements.begin(); it!=searchedElements.end(); it++)
            {
                unsigned int bulk_ele_idx = *it;
                ElementFullIter ele_2D = mesh->element(bulk_ele_idx);
                
                if (ele_2D->dim() == 2) { 
                    update_simplex(ele_2D, triangle_); // update triangle
                    
                    IntersectionAux<1,2> is(component_ele_idx, bulk_ele_idx, 0);
                    START_TIMER("Compute intersection");
                    ComputeIntersection<Simplex<1>, Simplex<2>> CI(abscissa_, triangle_);
                    CI.compute_final(is.points());
                    END_TIMER("Compute intersection");
                    
                    if(is.points().size() > 0) {
                        intersectionaux_storage12_.push_back(is);
                    }
                }
            }
            END_TIMER("Bounding box element iteration");
        }
    }

    END_TIMER("Element iteration");
}

// Declaration of specializations implemented in cpp:
template class InspectElementsAlgorithm<1>;
template class InspectElementsAlgorithm<2>;

} // END namespace
