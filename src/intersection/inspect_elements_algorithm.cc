/*
 *
 *
 *  Created on: 29.4.2016
 *      Author: pe
 */

#include <unordered_set>
#include <boost/functional/hash.hpp>

#include "inspect_elements_algorithm.hh"
#include "intersection_point_aux.hh"
#include "intersection_aux.hh"
#include "intersection_local.hh"
#include "compute_intersection.hh"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "mesh/ref_element.hh"
#include "mesh/bih_tree.hh"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"



template<unsigned int dimA, unsigned int dimB>
IntersectionAlgorithmBase<dimA,dimB>::IntersectionAlgorithmBase(Mesh* mesh)
: mesh(mesh)
{}

// template<unsigned int dimA, unsigned int dimB>
// template<unsigned int simplex_dim>
// void IntersectionAlgorithmBase<dimA,dimB>::update_simplex(const ElementAccessor<3>& element, Simplex< simplex_dim >& simplex)
// {
//     ASSERT(simplex_dim == element.dim());
//     arma::vec3 *field_of_points[simplex_dim+1];
//     for(unsigned int i=0; i < simplex_dim+1; i++)
//         field_of_points[i]= &(element.node(i)->point());
//     simplex.set_simplices(field_of_points);
// }


template<unsigned int dim>    
InspectElementsAlgorithm<dim>::InspectElementsAlgorithm(Mesh* input_mesh)
: IntersectionAlgorithmBase<dim,3>(input_mesh)
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
    END_TIMER("Intersection initialization");
}

template<unsigned int dim>
void InspectElementsAlgorithm<dim>::compute_bounding_boxes()
{
    START_TIMER("Compute bounding boxes");
    if(elements_bb.size() == 0){
        elements_bb.resize(mesh->n_elements());
        bool first_3d_element = true;
        for (auto elm : mesh->elements_range()) {

            elements_bb[elm.idx()] = elm.bounding_box();

                if (elm->dim() == 3){
                    if(first_3d_element){
                        first_3d_element = false;
                        mesh_3D_bb = elements_bb[elm.idx()];
                    }else{
                        mesh_3D_bb.expand(elements_bb[elm.idx()].min());
                        mesh_3D_bb.expand(elements_bb[elm.idx()].max());
                    }
                }
        }
    }
    END_TIMER("Compute bounding boxes");
}

template<unsigned int dim> 
bool InspectElementsAlgorithm<dim>::compute_initial_CI(const ElementAccessor<3> &comp_ele,
                                                       const ElementAccessor<3> &bulk_ele)
{
    unsigned int component_ele_idx = comp_ele.idx(),
                 bulk_ele_idx = bulk_ele.idx();
    
    IntersectionAux<dim,3> is(component_ele_idx, bulk_ele_idx);
    START_TIMER("Compute intersection");
    ComputeIntersection<dim,3> CI(comp_ele, bulk_ele, mesh);
    CI.init();
    CI.compute(is);
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
    //DebugOut() << "#########   ALGORITHM: compute_intersections   #########\n";
    
    init();
    
    START_TIMER("Element iteration");
    
    for (auto elm : mesh->elements_range()) {
        unsigned int component_ele_idx = elm.idx();
        
        if (elm->dim() == dim &&                                // is component element
            !closed_elements[component_ele_idx] &&                    // is not closed yet
            bih.ele_bounding_box(component_ele_idx).intersect(bih.tree_box()))    // its bounding box intersects 3D mesh bounding box
        {    
            std::vector<unsigned int> searchedElements;
            
            START_TIMER("BIHtree find");
            bih.find_bounding_box(bih.ele_bounding_box(component_ele_idx), searchedElements);
            END_TIMER("BIHtree find");

            START_TIMER("Bounding box element iteration");
            
            // Go through all element which bounding box intersects the component element bounding box
            for (std::vector<unsigned int>::iterator it = searchedElements.begin(); it!=searchedElements.end(); it++)
            {
                unsigned int bulk_ele_idx = *it;
                ElementAccessor<3> ele_3D = mesh->element_accessor( bulk_ele_idx );

                // if:
                // check 3D only
                // check with the last component element computed for the current 3D element
                // intersection has not been computed already
                if (ele_3D->dim() == 3 &&
                    (last_slave_for_3D_elements[bulk_ele_idx] != component_ele_idx &&
                     !intersection_exists(component_ele_idx,bulk_ele_idx) )
                ) {
                    // check that tetrahedron element is numbered correctly and is not degenerated
                    ASSERT_DBG(ele_3D.tetrahedron_jacobian() > 0).add_value(ele_3D.index(),"element index").error(
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
                    
                    bool found = compute_initial_CI(elm, ele_3D);

                    // keep the index of the current component element that is being investigated
                    unsigned int current_component_element_idx = component_ele_idx;
                    
                    if(found){
                        
                        prolongation_decide(elm, ele_3D, intersection_list_[component_ele_idx].back());
                        
                        START_TIMER("Prolongation algorithm");
                        do{
                            // flag is set false if the component element is not fully covered with tetrahedrons
                            bool element_covered = true;
                            
                            while(!bulk_queue_.empty()){
                                Prolongation pr = bulk_queue_.front();
                                //DebugOut().fmt("Bulk queue: ele_idx {}.\n",pr.elm_3D_idx);
                                
                                if( pr.elm_3D_idx == undefined_elm_idx_)
                                {
                                    //DebugOut().fmt("Open intersection component element: {}\n",current_component_element_idx);
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
                                //DebugOut().fmt("Component queue: ele_idx {}.\n",current_component_element_idx);
                                
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
    
    MessageOut().fmt("{}D-3D: number of intersections = {}\n", dim, n_intersections_);
    // DBG write which elements are closed
//     for (auto ele : mesh->elements_range()) {
//         DebugOut().fmt("Element[{}] closed: {}\n",ele.index(),(closed_elements[ele.index()] ? 1 : 0));
//     }
}
  
template<unsigned int dim>
void InspectElementsAlgorithm<dim>::compute_intersections_BIHtree(const BIHTree& bih)
{
    DebugOut() << "#########   ALGORITHM: compute_intersections_BIHtree   #########\n";
    
    init();
    
    START_TIMER("Element iteration");
    
    for (auto elm : mesh->elements_range()) {
        unsigned int component_ele_idx = elm.idx();
        
        if (elm.dim() == dim &&                                    // is component element
            bih.ele_bounding_box(component_ele_idx).intersect(bih.tree_box()))   // its bounding box intersects 3D mesh bounding box
        {   
            std::vector<unsigned int> searchedElements;
            
            START_TIMER("BIHtree find");
            bih.find_bounding_box(bih.ele_bounding_box(component_ele_idx), searchedElements);
            END_TIMER("BIHtree find");
            
            START_TIMER("Bounding box element iteration");
            
            // Go through all element which bounding box intersects the component element bounding box
            for (std::vector<unsigned int>::iterator it = searchedElements.begin(); it!=searchedElements.end(); it++)
            {
                unsigned int bulk_ele_idx = *it;
                ElementAccessor<3> ele_3D = mesh->element_accessor( bulk_ele_idx );
                
                if (ele_3D.dim() == 3
                ) {
                    // check that tetrahedron element is numbered correctly and is not degenerated
                    ASSERT_DBG(ele_3D.tetrahedron_jacobian() > 0).add_value(ele_3D.idx(),"element index").error(
                           "Tetrahedron element (%d) has wrong numbering or is degenerated (negative Jacobian).");
                    
                    IntersectionAux<dim,3> is(component_ele_idx, bulk_ele_idx);
                    START_TIMER("Compute intersection");
                    ComputeIntersection<dim,3> CI(elm, ele_3D, mesh);
                    CI.init();
                    CI.compute(is);
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
    //DebugOut() << "#########   ALGORITHM: compute_intersections_BB   #########\n";
    init();
    compute_bounding_boxes();
    
    START_TIMER("Element iteration");
    
    
    for (auto elm : mesh->elements_range()) {
        unsigned int component_ele_idx = elm.idx();
        
        if (elm.dim() == dim &&                                // is component element
            !closed_elements[component_ele_idx] &&                    // is not closed yet
            elements_bb[component_ele_idx].intersect(mesh_3D_bb))    // its bounding box intersects 3D mesh bounding box
        {    
            
            START_TIMER("Bounding box element iteration");
            
            // Go through all element which bounding box intersects the component element bounding box
            for (auto ele_3D : mesh->elements_range()) {
                unsigned int bulk_ele_idx = ele_3D.idx();

                // if:
                // check 3D only
                // check with the last component element computed for the current 3D element
                // check that the bounding boxes intersect
                // intersection has not been computed already
                if (ele_3D.dim() == 3 &&
                    (last_slave_for_3D_elements[bulk_ele_idx] != component_ele_idx &&
                     elements_bb[component_ele_idx].intersect(elements_bb[bulk_ele_idx]) &&
                     !intersection_exists(component_ele_idx,bulk_ele_idx) )
                ){
                    // check that tetrahedron element is numbered correctly and is not degenerated
                    ASSERT_DBG(ele_3D.tetrahedron_jacobian() > 0).add_value(ele_3D.index(),"element index").error(
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
                    
                    bool found = compute_initial_CI(elm, ele_3D);

                    // keep the index of the current component element that is being investigated
                    unsigned int current_component_element_idx = component_ele_idx;
                    
                    if(found){
                        //DebugOut().fmt("start component with elements {} {}\n",component_ele_idx, bulk_ele_idx);
                        
                        prolongation_decide(elm, ele_3D, intersection_list_[component_ele_idx].back());
                        
                        START_TIMER("Prolongation algorithm");
                        do{
                            // flag is set false if the component element is not fully covered with tetrahedrons
                            bool element_covered = true;
                            
                            while(!bulk_queue_.empty()){
                                Prolongation pr = bulk_queue_.front();
                                //DebugOut().fmt("Bulk queue: ele_idx {}.\n",pr.elm_3D_idx);
                                
                                if( pr.elm_3D_idx == undefined_elm_idx_)
                                {
                                    //DebugOut().fmt("Open intersection component element: {}\n",current_component_element_idx);
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
                                //DebugOut().fmt("Component queue: ele_idx {}.\n",current_component_element_idx);
                                
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
//     for (auto ele : mesh->elements_range()) {
//         DebugOut().fmt("Element[{}] closed: {}\n",ele.index(),(closed_elements[ele.index()] ? 1 : 0));
//     }
}



template<unsigned int dim>
template<unsigned int ele_dim>
std::vector< unsigned int > InspectElementsAlgorithm<dim>::get_element_neighbors(const ElementAccessor<3>& ele,
                                                                                 unsigned int ip_dim,
                                                                                 unsigned int ip_obj_idx)
{
    std::vector<Edge> edges;
    edges.reserve(ele_dim - ip_dim);  // reserve number of possible edges

    //DebugOut() << "dim " << ele_dim << ": ";
    switch (ip_dim)
    {
        // IP is at a node of tetrahedron; possible edges are from all connected sides (3)
        case 0: if(ele_dim == 1) {
                    edges.push_back(mesh->edge(ele->edge_idx(ip_obj_idx)));
                    break;
                }
                
                for(unsigned int j=0; j < RefElement<ele_dim>::n_sides_per_node; j++){
                    unsigned int local_edge = RefElement<ele_dim>::interact(Interaction<ele_dim-1,0>(ip_obj_idx))[j];
                    edges.push_back(mesh->edge(ele->edge_idx(local_edge)));
                }
                //DebugOut() << "prolong (node)\n";
                break;
        
        // IP is on a line of tetrahedron; possible edges are from all connected sides (2)
        case 1: if(ele_dim == 2) {
                    edges.push_back(mesh->edge(ele->edge_idx(ip_obj_idx)));
                    break;
                }
            
                ASSERT_DBG(ele_dim == 3);
                for(unsigned int j=0; j < RefElement<ele_dim>::n_sides_per_line; j++){
                    unsigned int local_edge = RefElement<ele_dim>::interact(Interaction<2,1>(ip_obj_idx))[j];
                    edges.push_back(mesh->edge(ele->edge_idx(local_edge)));
                }
                //DebugOut() << "prolong (edge)\n";
                break;
                
        // IP is on a side of tetrahedron; only possible edge is from the given side (1)
        case 2: ASSERT_DBG(ele_dim == 3);
                edges.push_back(mesh->edge(ele->edge_idx(ip_obj_idx)));
                //DebugOut() << "prolong (side)\n";
                break;
        default: ASSERT_DBG(0);
    }
    
    // get indices of neighboring bulk elements
    std::vector<unsigned int> elements_idx;
    elements_idx.reserve(2*edges.size());    // twice the number of edges
    for(Edge edg : edges)
    for(uint j=0; j < edg.n_sides();j++) {
        if ( edg.side(j)->element().idx() != ele.idx() )
            elements_idx.push_back(edg.side(j)->element().idx());
    }
    
    return elements_idx;
}


template<unsigned int dim>
unsigned int InspectElementsAlgorithm<dim>::create_prolongation(unsigned int bulk_ele_idx,
                                                                unsigned int component_ele_idx,
                                                                std::queue< Prolongation >& queue)
{
//     if(last_slave_for_3D_elements[bulk_ele_idx] == undefined_elm_idx_ ||
//         (last_slave_for_3D_elements[bulk_ele_idx] != component_ele_idx && !intersection_exists(component_ele_idx,bulk_ele_idx)))
//     {
    last_slave_for_3D_elements[bulk_ele_idx] = component_ele_idx;

    //DebugOut().fmt("prolongation: c {} in b {}\n",component_ele_idx,bulk_ele_idx);
    
    // prepare empty intersection object
    IntersectionAux<dim,3> il_other(component_ele_idx, bulk_ele_idx);
    intersection_list_[component_ele_idx].push_back(il_other);
    
    Prolongation pr = {component_ele_idx, bulk_ele_idx, (unsigned int)intersection_list_[component_ele_idx].size() - 1};
    queue.push(pr);
    
    return 1;
//     }
//     return 0;
}



template<unsigned int dim>
void InspectElementsAlgorithm<dim>::prolongation_decide(const ElementAccessor<3>& comp_ele,
                                                        const ElementAccessor<3>& bulk_ele,
                                                        IntersectionAux<dim,3>& is)
{
    //DebugOut() << "DECIDE\n";
    // number of IPs that are at vertices of component element (counter used for closing element)
    unsigned int n_ip_vertices = 0;
    
    if(dim == 2){
        // check whether all IPs lie in the face
        unsigned int sid = is.ips_in_face();
        if(sid < RefElement<3>::count<2>()){
            //NOTE: we are unable to combine cases A and B in 1d-2d
            // CASE A: compatible vb neighboring
//             DBGVAR(comp_ele.idx());
//             DBGVAR(bulk_ele.idx());
            if(comp_ele->n_neighs_vb() > 1){
                //set n_neighs_vb duplicities
                is.set_duplicities(comp_ele->n_neighs_vb());
//                 DBGCOUT(<< "vb neigh: N = " << is.duplicities() << "\n");
                //possibly copy intersection object
            }
            // CASE B: incompatible neighboring with lower dim element
            else {
                Edge edg = bulk_ele.side(sid)->edge();
                if(edg.n_sides() > 1){
                    //set n_sides duplicities
                    is.set_duplicities(edg.n_sides());
//                     DBGCOUT(<< "incomp. neigh: N = " << is.duplicities() << "\n");
                    //possibly copy intersection object
                }
            }
        }
    }
    
    for(const IntersectionPointAux<dim,3> &IP : is.points()) {
        
        // 1) prolong over component, if IP is at its boundary
        // 2) prolong over bulk, if IP is at its boundary
        // both cases are possible, so both prolongations might happen at once..
        
        if(IP.dim_A() < dim) { // if IP on the boundary of component element
            if(IP.dim_A() == 0) n_ip_vertices++;
            //DebugOut() << "on " << dim << "D boundary, dim = " << IP.dim_A() << "\n";
            
            // search for indices of neighboring component elements (including the current one)
            std::vector<unsigned int> comp_neighbors = get_element_neighbors<dim>(comp_ele,IP.dim_A(), IP.idx_A());
            
//             DBGCOUT( << comp_ele.idx() << "--" << bulk_ele.idx() << ":    ");
//             for(unsigned int& comp_neighbor_idx : comp_neighbors)
//                 cout << mesh->element_accessor(comp_neighbor_idx).idx() << "  ";
//             cout << "\n";
            
            unsigned int bulk_current = bulk_ele.idx();
            
            // add all component neighbors with current bulk element into component queue
            for(unsigned int& comp_neighbor_idx : comp_neighbors) {
                if(!intersection_exists(comp_neighbor_idx,bulk_current))
                    create_prolongation(bulk_current, comp_neighbor_idx, component_queue_);
            }
        }   
        
        if(IP.dim_B() < 3)
        {
            //DebugOut() << "on 3D boundary, dim = " << IP.dim_B() << "\n";
            
            // search for indices of neighboring bulk elements (including the current one)
            std::vector<unsigned int> bulk_neighbors = get_element_neighbors<3>(bulk_ele,IP.dim_B(),IP.idx_B());
            
//             DBGCOUT( << comp_ele.idx() << "--" << bulk_ele.idx() << ":    ");
//             for(unsigned int& bulk_neighbor_idx : bulk_neighbors)
//                 cout << mesh->element_accessor(bulk_neighbor_idx).idx() << "  ";
//             cout << "\n";
            
            unsigned int comp_current = comp_ele.idx();
            unsigned int n_prolongations = 0;
            // prolong over current comp element to other bulk elements (into bulk queue) (covering comp ele)
            for(unsigned int& bulk_neighbor_idx : bulk_neighbors)
            {
                if(last_slave_for_3D_elements[bulk_neighbor_idx] == undefined_elm_idx_ ||
                    (last_slave_for_3D_elements[bulk_neighbor_idx] != comp_current && 
                        !intersection_exists(comp_current,bulk_neighbor_idx)))
                    n_prolongations += create_prolongation(bulk_neighbor_idx,
                                                           comp_current,
                                                           bulk_queue_);
            }
            
            // if there are no sides of any edge that we can continue to prolongate over,
            // it means we are at the boundary and cannot prolongate further
            if(n_prolongations == 0)
            {
                Prolongation pr = {comp_ele.idx(), undefined_elm_idx_, undefined_elm_idx_};
                bulk_queue_.push(pr);
            }
        }
    }
    
    // close component element if it has all vertices inside bulk element
    if(n_ip_vertices == is.size()) closed_elements[comp_ele.idx()] = true;
}



template<>
void InspectElementsAlgorithm<2>::assert_same_intersection(unsigned int comp_ele_idx, unsigned int bulk_ele_idx)
{
    for(unsigned int i=0; i < intersection_list_[comp_ele_idx].size(); i++)
    {
        if(intersection_list_[comp_ele_idx][i].bulk_ele_idx() == bulk_ele_idx)
        {
            //DebugOut().fmt("intersection comp-bulk: {} {}\n", comp_ele_idx, bulk_ele_idx);
            ASSERT_DBG(0).add_value(bulk_ele_idx,"bulk_ele_idx").error("Want to add the same intersection!");
        }
    }
}


template<unsigned int dim>
void InspectElementsAlgorithm<dim>::prolongate(const InspectElementsAlgorithm< dim >::Prolongation& pr)
{
	ElementAccessor<3> elm = mesh->element_accessor( pr.component_elm_idx );
	ElementAccessor<3> ele_3D = mesh->element_accessor( pr.elm_3D_idx );
    
//     DebugOut().fmt("Prolongate: {} in {}.\n", elm.idx(), ele_3D.idx());

    //TODO: optimization: this might be called before and not every time 
    //(component element is not changing when emptying bulk queue)
//     this->update_simplex(elm, simplexA);
//     this->update_simplex(ele_3D, simplexB);

    IntersectionAux<dim,3> &is = intersection_list_[pr.component_elm_idx][pr.dictionary_idx];
    
    START_TIMER("Compute intersection");
    ComputeIntersection<dim,3> CI(elm, ele_3D, mesh);
    CI.init();
    CI.compute(is);
    END_TIMER("Compute intersection");
    
    last_slave_for_3D_elements[pr.elm_3D_idx] = pr.component_elm_idx;
    
    if(is.size() > 0){
//         for(unsigned int j=0; j < is.size(); j++) 
//             DebugOut() << is[j];
//         DebugOut().fmt("intersection of elements {} {} [{}--{}] size {}\n",
//                        elm.idx(), ele_3D.idx(),
//                        elm.region().label(), ele_3D.region().label(),
//                        is.size()
//                       );
        
        prolongation_decide(elm, ele_3D,is);
        n_intersections_++;
//         DBGVAR(n_intersections_);
    }
    else{
        // NOTE: we get here, when create_prolongation creates an empty intersection
        // - it can happen, that the CI will be empty
        // - currently, we remove it when storing the final intersection objects
        // - or we can erase it at this point
//         WarningOut() << "zero is: c " << elm.index() << " b " << ele_3D.index();
//         auto & v = intersection_list_[pr.component_elm_idx];
//         v.erase( v.next(v.being(),pr.dictionary_idx) );
    }
}




InspectElementsAlgorithm22::InspectElementsAlgorithm22(Mesh* input_mesh)
: IntersectionAlgorithmBase<2,2>(input_mesh)
{}


void InspectElementsAlgorithm22::compute_intersections(std::vector< std::vector<ILpair>>& intersection_map,
                                                       std::vector<IntersectionLocal<2,2>> &storage)
{
//     DebugOut() << "Intersections 2d-2d\n";
    ASSERT(storage.size() == 0);
    create_component_numbering();
    
    unsigned int ele_idx, eleA_idx, eleB_idx,
                 componentA_idx, componentB_idx,
                 temp_eleA_idx;
                 
    typedef std::pair<unsigned int, unsigned int> ipair;
    std::unordered_set<ipair, boost::hash<ipair>> computed_pairs;
    
    for (auto ele : mesh->elements_range()) {
    if (ele->dim() == 3)
    {
        ele_idx = ele.idx();
        // if there are not at least 2 2D elements intersecting 3D element; continue
        if(intersection_map[ele_idx].size() < 2) continue;
        
        const std::vector<ILpair> &local_map = intersection_map[ele_idx];
        
//         DebugOut() << print_var(local_map.size());
        for(unsigned int i=0; i < local_map.size(); i++)
        {
            //TODO: 1] compute all plucker coords at once
            //TODO: 2] pass plucker coords from 2d-3d
            
            eleA_idx = local_map[i].first;
            ElementAccessor<3> eleA = mesh->element_accessor( eleA_idx );
            if(eleA->dim() !=2 ) continue;  //skip other dimension intersection
            componentA_idx = component_idx_[eleA_idx];
            
//             IntersectionLocalBase * ilb = local_map[i].second;
//             DebugOut().fmt("2d-2d ILB: {} {} {}\n", ilb->bulk_ele_idx(), ilb->component_ele_idx(), componentA_idx);
            
            for(unsigned int j=i+1; j < local_map.size(); j++)
            {
                eleB_idx = local_map[j].first;
                componentB_idx = component_idx_[eleB_idx];

                if(componentA_idx == componentB_idx) continue;  //skip elements of the same component
                // this also skips the compatible connections (it is still a single component in this case)
                
                ElementAccessor<3> eleB = mesh->element_accessor( eleB_idx );
                if(eleB->dim() !=2 ) continue;  //skip other dimension intersection
                
                // set master -- slave order
                // do not overwrite the original eleA
                temp_eleA_idx = eleA_idx;
                ElementAccessor<3> temp_eleA = eleA;
                if (componentA_idx < componentB_idx){
                    std::swap(temp_eleA_idx, eleB_idx);
                    std::swap(temp_eleA, eleB);
                }
                
                //skip candidates already computed
                ipair ip = std::make_pair(temp_eleA_idx, eleB_idx);
                if(computed_pairs.count(ip) == 1){
//                     DBGCOUT(<< "skip: " << eleA_idx << " " << eleB_idx << "\n");
                    continue;
                }
                else{
                    compute_single_intersection(temp_eleA, eleB, storage);
                    computed_pairs.emplace(ip);
                }
                
//                 bool skip = false;
//                 for(unsigned int k=0; k<intersection_map[eleA_idx].size(); k++)
//                 {
//                     if(intersection_map[eleA_idx][k].first == eleB_idx) {
//                         skip = true;
//                         break;
//                     }
//                 }
//                 if(skip) continue;
                
//                 DebugOut().fmt("compute intersection 2d-2d: e_{} e_{} c_{} c_{}\n",
//                                eleA.index(), eleB.index(), componentA_idx, componentB_idx);
//                 DebugOut().fmt("compute intersection 2d-2d: e_{} e_{} c_{} c_{}\n",
//                                eleA.idx(), eleB.idx(), componentA_idx, componentB_idx);
                
//                 IntersectionAux<2,2> is;
//                 if (componentA_idx < componentB_idx)
//                     compute_single_intersection(eleA, eleB, storage);
//                 else
//                     compute_single_intersection(eleB, eleA, storage);
            }
        }
    }
    }
    MessageOut() << "2D-2D: number of intersections = " << storage.size() << "\n";
}


void InspectElementsAlgorithm22::compute_single_intersection(const ElementAccessor<3>& eleA,
                                                             const ElementAccessor<3>& eleB,
                                                             std::vector<IntersectionLocal<2,2>> &storage)
{
    ASSERT_DBG(eleA.dim() == 2);
    ASSERT_DBG(eleB.dim() == 2);
    ASSERT_DBG(eleA.idx() != eleB.idx());
    
    IntersectionAux<2,2> is(eleA.idx(), eleB.idx());
    
    ComputeIntersection<2,2> CI(eleA, eleB, mesh);
    CI.init();
    unsigned int n_local_intersection = CI.compute(is);
    
    // do not store point intersections
    if(n_local_intersection > 1){
        storage.push_back(IntersectionLocal<2,2>(is));
    }

}

void InspectElementsAlgorithm22::create_component_numbering()
{
    component_idx_.resize(mesh->n_elements(),unset_comp);
    component_counter_ = 0;
    
    // prolongation queue in the component mesh.
    std::queue<unsigned int> queue;

    for (auto ele : mesh->elements_range()) {
        if (ele->dim() == 2 &&
            component_idx_[ele.idx()] == (unsigned int)-1)
        {
            // start component
            queue.push(ele.idx());
            
            while(!queue.empty()){
                unsigned int ele_idx = queue.front();
                queue.pop();
                const ElementAccessor<3>& elm = mesh->element_accessor( ele_idx );
                for(unsigned int sid=0; sid < elm->n_sides(); sid++) {
                    Edge edg = elm.side(sid)->edge();

                    for(uint j=0; j < edg.n_sides();j++) {
                        uint neigh_idx = edg.side(j)->element().idx();
                        if (component_idx_[neigh_idx] == (unsigned int)-1) {
                            component_idx_[neigh_idx] = component_counter_;
                            queue.push(neigh_idx);
                        }
                    }
                }
            }
            component_counter_++;
        }
    }
    
    MessageOut() << "2D-2D: number of components = " << component_counter_ << "\n";
    
//     DBGCOUT(<< "Component numbering: \n");
//     for (auto ele : mesh->elements_range()) {
//         if (ele->dim() == 2){
//             cout << "2d ele " << ele.index() << ":  " << component_idx_[ele.index()] << endl;
//         }
//     }
}





InspectElementsAlgorithm12::InspectElementsAlgorithm12(Mesh* input_mesh)
: IntersectionAlgorithmBase<1,2>(input_mesh)
{}


void InspectElementsAlgorithm12::compute_intersections_3(std::vector< std::vector<ILpair>>& intersection_map,
                                                         std::vector<IntersectionLocal<1,2>> &storage)
{
    //DebugOut() << "Intersections 1d-2d\n";
    intersectionaux_storage12_.clear();
    ASSERT(storage.size() == 0);
    
    for (auto ele : mesh->elements_range()) {
    if (ele->dim() == 3)
    {
        unsigned int ele_idx = ele.idx();
        // if there are not at least 2 elements intersecting 3D element; continue
        if(intersection_map[ele_idx].size() < 2) continue;
        
        const std::vector<ILpair> &local_map = intersection_map[ele_idx];
        
        //DebugOut() << "more than 2 intersections in tetrahedron found\n";
        for(unsigned int i=0; i < local_map.size(); i++)
        {
            //TODO: 1] compute all plucker coords at once
            //TODO: 2] pass plucker coords from 1d-3d
            
            unsigned int eleA_idx = local_map[i].first;
            ElementAccessor<3> eleA = mesh->element_accessor( eleA_idx );
            
            if(eleA.dim() !=1 ) continue;  //skip other dimension intersection
            
            for(unsigned int j=0; j < local_map.size(); j++)
            {
                unsigned int eleB_idx = local_map[j].first;
                ElementAccessor<3> eleB = mesh->element_accessor( local_map[j].first );
                if(eleB.dim() !=2 ) continue;  //skip other dimension intersection
                
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
                
                //DebugOut().fmt("compute intersection 1d-2d: {} {}\n",eleA.index(), eleB.index());
//                 compute_single_intersection(eleA,
//                                             eleB);
                
                IntersectionAux<1,2> is(eleA_idx, eleB_idx);
                
                ComputeIntersection<1,2> CI(eleA, eleB, mesh);
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
                    
                    //DebugOut().fmt("1D-2D intersection [{} - {}]:\n",is.component_ele_idx(), is.bulk_ele_idx());
                }
            }
        }
    }
    }
    
    MessageOut() << "1D-2D [3]: number of intersections = " << storage.size() << "\n";
    // just dbg output
//     for(IntersectionLocal<1,2> &is : storage)
//     {
//         DebugOut().fmt("1D-2D intersection [{} - {}]:\n",is.component_ele_idx(), is.bulk_ele_idx());
//         for(const IntersectionPoint<1,2>& ip : is.points()) {
//             DebugOut() << ip;
//             auto p = ip.coords(mesh->element(is.component_ele_idx()));
//             DebugOut() << "[" << p[0] << " " << p[1] << " " << p[2] << "]\n";
//         }
//     }
}

// void InspectElementsAlgorithm12::compute_single_intersection(const ElementAccessor<3>& eleA,
//                                                              const ElementAccessor<3>& eleB)
// {
//     ASSERT_DBG(eleA.dim() == 1);
//     ASSERT_DBG(eleB.dim() == 2);
//     
//     this->update_simplex(eleA, simplexA);
//     this->update_simplex(eleB, simplexB);
//     
//     IntersectionAux<1,2> is(eleA.index(), eleB.index(), 0);
// //     std::vector<unsigned int> prolongation_table;
//     
//     ComputeIntersection< Simplex<1>, Simplex<2>> CI(simplexA, simplexB);
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
    //DebugOut() << "Intersections 1d-2d (2-bihtree)\n";
    intersectionaux_storage12_.clear();
    START_TIMER("Element iteration");
    
    for (auto elm : mesh->elements_range()) {
        unsigned int component_ele_idx = elm.idx();
        
        if (elm.dim() == 1)                                    // is component element
            //&& elements_bb[component_ele_idx].intersect(mesh_3D_bb))   // its bounding box intersects 3D mesh bounding box
        {   
            std::vector<unsigned int> searchedElements;
            
            START_TIMER("BIHtree find");
            bih.find_bounding_box(bih.ele_bounding_box(component_ele_idx), searchedElements);
            END_TIMER("BIHtree find");
            
            START_TIMER("Bounding box element iteration");
            
            // Go through all element which bounding box intersects the component element bounding box
            for (std::vector<unsigned int>::iterator it = searchedElements.begin(); it!=searchedElements.end(); it++)
            {
                unsigned int bulk_ele_idx = *it;
                ElementAccessor<3> ele_2D = mesh->element_accessor( bulk_ele_idx );
                
                if (ele_2D.dim() == 2) {
                    
                    IntersectionAux<1,2> is(component_ele_idx, bulk_ele_idx);
                    START_TIMER("Compute intersection");
                    ComputeIntersection<1,2> CI(elm, ele_2D, mesh);
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

void InspectElementsAlgorithm12::compute_intersections_1(const BIHTree& bih)
{
    //DebugOut() << "Intersections 1d-2d (2-bihtree) in 2D plane.\n";
    intersectionaux_storage12_.clear();
    START_TIMER("Element iteration");
    
    for (auto elm : mesh->elements_range()) {
        unsigned int component_ele_idx = elm.idx();
        
        if (elm->dim() == 1)                                    // is component element
            //&& elements_bb[component_ele_idx].intersect(mesh_3D_bb))   // its bounding box intersects 3D mesh bounding box
        {   
            std::vector<unsigned int> candidate_list;
            bih.find_bounding_box(bih.ele_bounding_box(component_ele_idx), candidate_list);
            
            START_TIMER("Bounding box element iteration");
            
            // Go through all element which bounding box intersects the component element bounding box
            for(unsigned int bulk_ele_idx : candidate_list) {
            	ElementAccessor<3> ele_2D = mesh->element_accessor(bulk_ele_idx);
                
                if (ele_2D->dim() == 2) { 
                    
                    IntersectionAux<1,2> is(component_ele_idx, bulk_ele_idx);
                    START_TIMER("Compute intersection");
                    ComputeIntersection<1,2> CI(elm, ele_2D, mesh);
                    CI.compute_final_in_plane(is.points());
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
    MessageOut() << "1D-2D [1]: number of intersections = " << intersectionaux_storage12_.size() << "\n";
}

// Declaration of specializations implemented in cpp:
template class IntersectionAlgorithmBase<1,3>;
template class IntersectionAlgorithmBase<2,3>;
template class IntersectionAlgorithmBase<1,2>;
template class IntersectionAlgorithmBase<2,2>;

template class InspectElementsAlgorithm<1>;
template class InspectElementsAlgorithm<2>;


