/*
 * inspectelements.cpp
 *
 *  Created on: 13.4.2014
 *      Author: viktor, pe, jb
 */

#include "inspectelements.h"
#include "intersectionpoint.h"
#include "intersectionaux.h"
#include "intersection_local.h"
#include "computeintersection.h"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "mesh/ref_element.hh"
#include "mesh/bih_tree.hh"

namespace computeintersection {

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
    END_TIMER("Intersection initialization");
}

template<unsigned int dim>
template<unsigned int simplex_dim>
void InspectElementsAlgorithm<dim>::update_simplex(const ElementFullIter& element, Simplex< simplex_dim >& simplex)
{
    arma::vec3 *field_of_points[simplex_dim+1];
    for(unsigned int i=0; i < simplex_dim+1; i++)
        field_of_points[i]= &(element->node[i]->point());
    simplex.set_simplices(field_of_points);
}

 
template<unsigned int dim> 
bool InspectElementsAlgorithm<dim>::compute_initial_CI(unsigned int component_ele_idx, unsigned int bulk_ele_idx,
                                                       std::vector<unsigned int> &prolongation_table)
{
    IntersectionAux<dim,3> is(component_ele_idx, bulk_ele_idx);
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
    bool found = false;

    for(unsigned int i = 0; i < intersection_list_[component_ele_idx].size();i++){

        if(intersection_list_[component_ele_idx][i].bulk_ele_idx() == bulk_ele_idx){
            found = true;
            break;
        }
    }
    return found;
}


template<unsigned int dim>
void InspectElementsAlgorithm<dim>::compute_intersections()
{
    init();
    START_TIMER("BIHTree");
    BIHTree bt(mesh, 20);
    END_TIMER("BIHTree");
    
    START_TIMER("Element iteration");
    
    FOR_ELEMENTS(mesh, elm) {
        unsigned int component_ele_idx = elm->index();
        
        if (elm->dim() == dim &&                                // is component element
            !closed_elements[component_ele_idx] &&                    // is not closed yet
            elements_bb[component_ele_idx].intersect(mesh_3D_bb))    // its bounding box intersects 3D mesh bounding box
        {    
            std::vector<unsigned int> searchedElements;
            bt.find_bounding_box(elements_bb[component_ele_idx], searchedElements);

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
                    (last_slave_for_3D_elements[bulk_ele_idx] != (int)(component_ele_idx) &&
                     !intersection_exists(component_ele_idx,bulk_ele_idx) )
                ) {
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
                    bool found = compute_initial_CI(component_ele_idx, bulk_ele_idx, prolongation_table);

                    // keep the index of the current component element that is being investigated
                    unsigned int current_component_element_idx = component_ele_idx;
                    
                    if(found){
                        DBGMSG("start component with elements %d %d\n",component_ele_idx, bulk_ele_idx);
                        
                        prolongation_decide(elm, ele_3D, intersection_list_[component_ele_idx].back(), prolongation_table);
                        
                        START_TIMER("Prolongation algorithm");
                        do{
                            // flag is set false if the component element is not fully covered with tetrahedrons
                            bool element_covered = true;
                            
                            while(!bulk_queue_.empty()){
                                Prolongation pr = bulk_queue_.front();
                                DBGMSG("Bulk queue: ele_idx %d.\n",pr.elm_3D_idx);
                                
                                if( pr.elm_3D_idx == undefined_elm_idx_)
                                {
                                    DBGMSG("Open intersection component element: %d\n",current_component_element_idx);
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
                                DBGMSG("Component queue: ele_idx %d.\n",current_component_element_idx);
                                
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
    FOR_ELEMENTS(mesh, ele) {
        DBGMSG("Element[%3d] closed: %d\n",ele.index(),(closed_elements[ele.index()] ? 1 : 0));
    }
}
  

template<unsigned int dim>
std::vector< unsigned int > InspectElementsAlgorithm<dim>::get_bulk_element_edges(const ElementFullIter& bulk_ele,
                                                                                  const IntersectionPointAux< dim, 3  >& IP,
                                                                                  const bool &include_current_bulk_ele
                                                                                 )
{
    std::vector<Edge*> edges;
    edges.reserve(3 - IP.dim_B());  // reserve number of possible edges

    switch (IP.dim_B())
    {
        // IP is at a node of tetrahedron; possible edges are from all connected sides (3)
        case 0: for(unsigned int j=0; j < RefElement<3>::n_sides_per_node; j++)
                    edges.push_back(&(mesh->edges[bulk_ele->edge_idx_[RefElement<3>::interact<2,0>(IP.idx_B())[j]]]));
                DBGMSG("3d prolong (node)\n");
                break;
        
        // IP is on a line of tetrahedron; possible edges are from all connected sides (2)
        case 1: for(unsigned int j=0; j < RefElement<3>::n_sides_per_line; j++)
                    edges.push_back(&(mesh->edges[bulk_ele->edge_idx_[RefElement<3>::interact<2,1>(IP.idx_B())[j]]]));
                DBGMSG("3d prolong (edge)\n");
                break;
                
        // IP is on a side of tetrahedron; only possible edge is from the given side (1)
        case 2: edges.push_back(&(mesh->edges[bulk_ele->edge_idx_[IP.idx_B()]]));
                DBGMSG("3d prolong (side)\n");
                break;
        default: ASSERT_LESS(IP.dim_B(),3);
    }
    
    // get indices of neighboring bulk elements
    std::vector<unsigned int> bulk_elements_idx;
    bulk_elements_idx.reserve(2*(3-IP.dim_B()));    // twice the number of edges
    for(Edge* edg : edges)
    for(int j=0; j < edg->n_sides;j++) {
        if (edg->side(j)->element() != bulk_ele)
            bulk_elements_idx.push_back(edg->side(j)->element()->index());
    }
    
    // possibly include the current bulk element
    if(include_current_bulk_ele)
        bulk_elements_idx.push_back(bulk_ele->index());
    
    return bulk_elements_idx;
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
        
            DBGMSG("3d prolong %d in %d\n",component_ele_idx,bulk_neighbor_idx);
            
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



template<>
void InspectElementsAlgorithm<1>::prolongation_decide(const ElementFullIter& comp_ele,
                                                      const ElementFullIter& bulk_ele,
                                                      const IntersectionAux<1,3> &is,
                                                      const std::vector<unsigned int> &prolongation_table)
{
    DBGMSG("DECIDE\n");
    // number of IPs that are at vertices of component element (counter used for closing element)
    unsigned int n_ip_vertices = 0;
    
    for(const IntersectionPointAux<1,3> &IP : is.points()) {
        if(IP.dim_A() == 0) { // if IP is the end of the 1D element
            n_ip_vertices++;
            DBGMSG("1D end\n");
            
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
                            DBGMSG("1d prolong %d in %d\n", component_neighbor_idx, bulk_ele->index());
                            
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
                std::vector<unsigned int> bulk_neighbors = get_bulk_element_edges(bulk_ele,IP, true);
                
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
            std::vector<unsigned int> bulk_neighbors = get_bulk_element_edges(bulk_ele,IP,false);
            
            // if edge has only one side, it means it is on the boundary and we cannot prolongate
            if(bulk_neighbors.empty())
            {
                Prolongation pr = {comp_ele->index(), undefined_elm_idx_, undefined_elm_idx_};
                bulk_queue_.push(pr);
                continue;
            }

            unsigned int n_prolongations = create_prolongations_over_bulk_element_edges(bulk_neighbors,comp_ele->index());
            
            DBGMSG("cover: %d %d\n", is.size(), n_prolongations);
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
void InspectElementsAlgorithm<2>::prolongation_decide(const ElementFullIter& comp_ele,
                                                      const ElementFullIter& bulk_ele,
                                                      const IntersectionAux<2,3> &is,
                                                      const std::vector<unsigned int> &prolongation_table
                                                     )
{
    for(unsigned int i = 0; i < prolongation_table.size();i++){

        unsigned int side;
        bool is_triangle_side = true;

        if(prolongation_table[i] >= 4){
            side = prolongation_table[i] - 4;
        }else{
            side = prolongation_table[i];
            is_triangle_side = false;
        }

        DBGMSG("prolongation table: %d %d\n", side, is_triangle_side);

        if(is_triangle_side){
            // prolongation through the triangle side

            SideIter elm_side = comp_ele->side(side);
            Edge *edg = elm_side->edge();

            for(int j=0; j < edg->n_sides;j++) {
                SideIter other_side=edg->side(j);
                if (other_side != elm_side) {
                    unsigned int sousedni_element = other_side->element()->index(); // 2D element

                    DBGMSG("2d sousedni_element %d\n", sousedni_element);

                        if(!intersection_exists(sousedni_element,bulk_ele->index())){

                            DBGMSG("2d prolong\n");
                            // Vytvoření průniku bez potřeby počítání
                            IntersectionAux<2,3> il_other(sousedni_element, bulk_ele->index());
                            intersection_list_[sousedni_element].push_back(il_other);

                            Prolongation pr = {sousedni_element, bulk_ele->index(), (unsigned int)intersection_list_[sousedni_element].size() - 1};
                            component_queue_.push(pr);
                        }
                }
            }

        }else{
            // prolongation through the tetrahedron side
            
            SideIter elm_side = bulk_ele->side(side);
            Edge *edg = elm_side->edge();

            for(int j=0; j < edg->n_sides;j++) {
                SideIter other_side=edg->side(j);
                if (other_side != elm_side) {
                    //

                    unsigned int sousedni_element = other_side->element()->index();

                    DBGMSG("3d sousedni_element %d\n", sousedni_element);

                    // TODO:
                    // - rename it, describe it and test that it is really useful !!
                    // how it probably works: 
                    // - if index undefined then no intersection has been computed for the 3D element
                    // - if not undefined check the index of actual 2D element 
                    //   (most probable case, that we are looking along 2D element back to 3D element, which we have just computed)
                    // - if it is another 2D element, then go through all found intersections of the 3D element and test it..
                    
                    if(last_slave_for_3D_elements[sousedni_element] == undefined_elm_idx_ || 
                        (last_slave_for_3D_elements[sousedni_element] != comp_ele->index() && !intersection_exists(comp_ele->index(),sousedni_element))){
                        
//                         last_slave_for_3D_elements[sousedni_element] = comp_ele->index();

                        DBGMSG("3d prolong\n");
                        
                        // Vytvoření průniku bez potřeby počítání
                        IntersectionAux<2,3> il_other(comp_ele->index(), sousedni_element);
                        intersection_list_[comp_ele->index()].push_back(il_other);

                        Prolongation pr = {comp_ele->index(), sousedni_element, (unsigned int)intersection_list_[comp_ele->index()].size() - 1};
                        bulk_queue_.push(pr);

                    }
                }
            }

        }

    }
}

template<unsigned int dim>
void InspectElementsAlgorithm<dim>::prolongate(const InspectElementsAlgorithm< dim >::Prolongation& pr)
{
    ElementFullIter elm = mesh->element(pr.component_elm_idx);
    ElementFullIter ele_3D = mesh->element(pr.elm_3D_idx);
    
    DBGMSG("Prolongate %dD: %d in %d.\n", dim, pr.component_elm_idx, pr.elm_3D_idx);

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
//             cout << is[j];
        
        prolongation_decide(elm, ele_3D,is,prolongation_table);
        n_intersections_++;
    }
}


 
 
InspectElements::InspectElements(Mesh* mesh)
: mesh(mesh), algorithm13_(mesh), algorithm23_(mesh), algorithm22_(mesh)
{}

InspectElements::~InspectElements()
{}

 
double InspectElements::measure_13() 
{
    double subtotal = 0.0;

    for(unsigned int i = 0; i < intersection_storage13_.size(); i++){
        ElementFullIter ele = mesh->element(intersection_storage13_[i].component_ele_idx());            
        double t1d_length = ele->measure();
        double local_length = intersection_storage13_[i].compute_measure();
        
        if(intersection_storage13_[i].size() == 2)
        {
        arma::vec3 from = intersection_storage13_[i][0].coords(ele);
        arma::vec3 to = intersection_storage13_[i][1].coords(ele);
        DBGMSG("sublength from [%f %f %f] to [%f %f %f] = %f\n",
               from[0], from[1], from[2], 
               to[0], to[1], to[2],
               local_length*t1d_length);
        }
        subtotal += local_length*t1d_length;
    }
    return subtotal;
}


double InspectElements::measure_23() 
{
    double subtotal = 0.0;

    for(unsigned int i = 0; i < intersection_storage23_.size(); i++){
            double t2dArea = mesh->element(intersection_storage23_[i].component_ele_idx())->measure();
            double localArea = intersection_storage23_[i].compute_measure();
            subtotal += 2*localArea*t2dArea;
        }
    return subtotal;
} 
 

template<unsigned int dim>
void InspectElements::compute_intersections(InspectElementsAlgorithm< dim >& iea,
                                            std::vector< IntersectionLocal<dim,3>>& storage)
{
    START_TIMER("Intersection algorithm");
    iea.compute_intersections();
    END_TIMER("Intersection algorithm");
    
    START_TIMER("Intersection into storage");
    storage.reserve(iea.n_intersections_);
    
    FOR_ELEMENTS(mesh, elm) {
        unsigned int idx = elm->index(); 
        unsigned int bulk_idx;
        
        if(elm->dim() == dim)
        {
                intersection_map_[idx].resize(iea.intersection_list_[idx].size());
                for(unsigned int j = 0; j < iea.intersection_list_[idx].size(); j++){
                    bulk_idx = iea.intersection_list_[idx][j].bulk_ele_idx();
//                     DBGMSG("cidx %d bidx %d: n=%d\n",idx, bulk_idx, iea_13.intersection_list_[idx][j].size());
                    storage.push_back(IntersectionLocal<dim,3>(iea.intersection_list_[idx][j]));
                    
                    // create map for component element
                    intersection_map_[idx][j] = std::make_pair(
                                                    bulk_idx, 
                                                    &(storage.back())
                                                );
//                  // write down intersections
//                     IntersectionLocal<1,3>* il13 = 
//                         static_cast<IntersectionLocal<1,3>*> (intersection_map_[idx][j].second);
//                     cout << &(intersection_storage13_.back()) << "  " << il13 << *il13;
//                     for(IntersectionPoint<1,3> &ip : il13->points())
//                         ip.coords(mesh->element(idx)).print(cout);

                    // create map for bulk element
                    intersection_map_[bulk_idx].push_back(
                                                std::make_pair(
                                                    idx, 
                                                    &(intersection_storage13_.back())
                                                ));
                }
        }
    }
    END_TIMER("Intersection into storage");
}

void InspectElements::compute_intersections_22(vector< IntersectionLocal< 2, 2 > >& storage)
{
    START_TIMER("Intersection algorithm");
    algorithm22_.compute_intersections(intersection_map_);
    END_TIMER("Intersection algorithm");
    
    START_TIMER("Intersection into storage");
    storage.reserve(algorithm22_.intersectionaux_storage22_.size());
    
    for(IntersectionAux<2,2> &is : algorithm22_.intersectionaux_storage22_) {
        unsigned int triaA_idx = is.component_ele_idx();
        unsigned int triaB_idx = is.bulk_ele_idx();

        storage.push_back(IntersectionLocal<2,2>(is));
        intersection_map_[triaA_idx].push_back(std::make_pair(
                                                    triaB_idx,
                                                    &(storage.back())
                                                ));
        intersection_map_[triaB_idx].push_back(std::make_pair(
                                                    triaA_idx,
                                                    &(storage.back())
                                                ));
    }
    END_TIMER("Intersection into storage");
}

void InspectElements::compute_intersections(computeintersection::IntersectionType d)
{
    intersection_map_.resize(mesh->n_elements());
    
    if(d & IntersectionType::d13){
        START_TIMER("Intersections 1D-3D");
        compute_intersections<1>(algorithm13_,intersection_storage13_);
        END_TIMER("Intersections 1D-3D");
    }
    
    if(d & IntersectionType::d23 | IntersectionType::d22){
        START_TIMER("Intersections 2D-3D");
        compute_intersections<2>(algorithm23_,intersection_storage23_);
        END_TIMER("Intersections 2D-3D");
    }
    
     if(d & IntersectionType::d22){
        START_TIMER("Intersections 2D-2D");
        compute_intersections_22(intersection_storage22_);
        END_TIMER("Intersections 2D-2D");
    }
}

 
void InspectElements::print_mesh_to_file_13(string name)
{
        string t_name = name;

        unsigned int number_of_intersection_points = 0;
        unsigned int number_of_nodes = mesh->n_nodes();

        for(unsigned int j = 0; j < intersection_storage13_.size();j++){
                number_of_intersection_points += intersection_storage13_[j].size();
        }

        FILE * file;
        file = fopen((t_name.append(".msh")).c_str(),"w");

        fprintf(file, "$MeshFormat\n");
        fprintf(file, "2.2 0 8\n");
        fprintf(file, "$EndMeshFormat\n");
        fprintf(file, "$Nodes\n");
        fprintf(file, "%d\n", (mesh->n_nodes() + number_of_intersection_points));

        FOR_NODES(mesh, nod){
            arma::vec3 _nod = nod->point();
            fprintf(file,"%d %f %f %f\n", nod.id(), _nod[0], _nod[1], _nod[2]);
        }

        for(unsigned int j = 0; j < intersection_storage13_.size();j++){
            IntersectionLocal<1,3> il = intersection_storage13_[j];
            ElementFullIter el1D = mesh->element(il.component_ele_idx());
//             ElementFullIter el3D = mesh->element(il.bulk_ele_idx());
            
            for(unsigned int k = 0; k < il.size();k++){
                number_of_nodes++;
                IntersectionPoint<1,3> IP13 = il[k];
                arma::vec3 global = IP13.coords(el1D);
                
//                 if(i == 0){
//                     _global = (IP13.local_bcoords_A())[0] * el1D->node[0]->point()
//                                                        +(IP13.local_bcoords_A())[1] * el1D->node[1]->point();
//                 }else{
//                     _global = (IP13.local_bcoords_B())[0] * el3D->node[0]->point()
//                                                        +(IP13.local_bcoords_B())[1] * el3D->node[1]->point()
//                                                        +(IP13.local_bcoords_B())[2] * el3D->node[2]->point()
//                                                        +(IP13.local_bcoords_B())[3] * el3D->node[3]->point();
//                 }

                fprintf(file,"%d %f %f %f\n", number_of_nodes, global[0], global[1], global[2]);
            }
        }

        fprintf(file,"$EndNodes\n");
        fprintf(file,"$Elements\n");
        fprintf(file,"%d\n", ((unsigned int)intersection_storage13_.size() + mesh->n_elements()) );

        FOR_ELEMENTS(mesh, elee){
            if(elee->dim() == 3){
                int id1 = mesh->node_vector.index(elee->node[0]) + 1;
                int id2 = mesh->node_vector.index(elee->node[1]) + 1;
                int id3 = mesh->node_vector.index(elee->node[2]) + 1;
                int id4 = mesh->node_vector.index(elee->node[3]) + 1;

                fprintf(file,"%d 4 2 %d %d %d %d %d %d\n", elee.id(), elee->region().id(), elee->pid, id1, id2, id3, id4);
            }else if(elee->dim() == 2){
                int id1 = mesh->node_vector.index(elee->node[0]) + 1;
                int id2 = mesh->node_vector.index(elee->node[1]) + 1;
                int id3 = mesh->node_vector.index(elee->node[2]) + 1;
                fprintf(file,"%d 2 2 %d %d %d %d %d\n", elee.id(), elee->region().id(), elee->pid, id1, id2, id3);

            }else if(elee->dim() == 1){
                int id1 = mesh->node_vector.index(elee->node[0]) + 1;
                int id2 = mesh->node_vector.index(elee->node[1]) + 1;
                fprintf(file,"%d 1 2 %d %d %d %d\n",elee.id(), elee->region().id(), elee->pid, id1, id2);
            }
        }

        unsigned int number_of_elements = mesh->n_elements();
        unsigned int nodes = mesh->n_nodes();

        for(unsigned int j = 0; j < intersection_storage13_.size();j++){
            IntersectionLocal<1,3> il = intersection_storage13_[j];
            number_of_elements++;
            nodes++;
            if(il.size() == 1){
                fprintf(file,"%d 1 2 1001 0 %d %d\n", number_of_elements, nodes, nodes);
            }else if(il.size() == 2){
                fprintf(file,"%d 1 2 1001 0 %d %d\n", number_of_elements, nodes, nodes+1);
                nodes++;
            }
        }

        fprintf(file,"$EndElements\n");
        fclose(file);
}

void InspectElements::print_mesh_to_file_23(string name)
{
    //for(unsigned int i = 0; i < 2;i++){
        string t_name = name;

        unsigned int number_of_intersection_points = 0;
        unsigned int number_of_nodes = mesh->n_nodes();

        for(unsigned int j = 0; j < intersection_storage23_.size();j++){
            number_of_intersection_points += intersection_storage23_[j].size();
        }

        FILE * file;
        file = fopen((t_name.append(".msh")).c_str(),"w");
        
        fprintf(file, "$MeshFormat\n");
        fprintf(file, "2.2 0 8\n");
        fprintf(file, "$EndMeshFormat\n");
        fprintf(file, "$Nodes\n");
        fprintf(file, "%d\n", (mesh->n_nodes() + number_of_intersection_points));

        FOR_NODES(mesh, nod){
            arma::vec3 _nod = nod->point();
            fprintf(file,"%d %f %f %f\n", nod.id(), _nod[0], _nod[1], _nod[2]);
        }

        for(unsigned int j = 0; j < intersection_storage23_.size();j++){
            
            IntersectionLocal<2,3> il = intersection_storage23_[j];
            ElementFullIter el2D = mesh->element(il.component_ele_idx());
//             ElementFullIter el3D = mesh->element(il.bulk_ele_idx());
            
            for(unsigned int k = 0; k < intersection_storage23_[j].size();k++){

                    number_of_nodes++;
                    IntersectionPoint<2,3> IP23 = il[k];
                    arma::vec3 global = IP23.coords(el2D);
//                     if(i == 0){
//                         _global = (IP23.local_bcoords_A())[0] * el2D->node[0]->point()
//                                                            +(IP23.local_bcoords_A())[1] * el2D->node[1]->point()
//                                                            +(IP23.local_bcoords_A())[2] * el2D->node[2]->point();
//                     }else{
//                         _global = (IP23.local_bcoords_B())[0] * el3D->node[0]->point()
//                                                            +(IP23.local_bcoords_B())[1] * el3D->node[1]->point()
//                                                            +(IP23.local_bcoords_B())[2] * el3D->node[2]->point()
//                                                            +(IP23.local_bcoords_B())[3] * el3D->node[3]->point();
//                     }
                    fprintf(file,"%d %f %f %f\n", number_of_nodes, global[0], global[1], global[2]);
            }
        }

        fprintf(file,"$EndNodes\n");
        fprintf(file,"$Elements\n");
        fprintf(file,"%d\n", (number_of_intersection_points + mesh->n_elements()) );

        FOR_ELEMENTS(mesh, elee){
            if(elee->dim() == 3){
                int id1 = mesh->node_vector.index(elee->node[0]) + 1;
                int id2 = mesh->node_vector.index(elee->node[1]) + 1;
                int id3 = mesh->node_vector.index(elee->node[2]) + 1;
                int id4 = mesh->node_vector.index(elee->node[3]) + 1;

                fprintf(file,"%d 4 2 %d %d %d %d %d %d\n", elee.id(), elee->region().id(), elee->pid, id1, id2, id3, id4);
            }else if(elee->dim() == 2){
                int id1 = mesh->node_vector.index(elee->node[0]) + 1;
                int id2 = mesh->node_vector.index(elee->node[1]) + 1;
                int id3 = mesh->node_vector.index(elee->node[2]) + 1;
                fprintf(file,"%d 2 2 %d %d %d %d %d\n", elee.id(), elee->region().id(), elee->pid, id1, id2, id3);

            }else{
                int id1 = mesh->node_vector.index(elee->node[0]) + 1;
                int id2 = mesh->node_vector.index(elee->node[1]) + 1;
                fprintf(file,"%d 1 2 %d %d %d %d\n",elee.id(), elee->region().id(), elee->pid, id1, id2);
            }
        }

        unsigned int number_of_elements = mesh->n_elements();
        unsigned int last = 0;
        unsigned int nodes = mesh->n_nodes();

        for(unsigned int j = 0; j < intersection_storage23_.size();j++){
            IntersectionLocal<2,3> il = intersection_storage23_[j];
            
                for(unsigned int k = 0; k < il.size();k++){

                    number_of_elements++;
                    nodes++;
                    if(k == 0){
                        last = nodes;
                    }

                    if((k+1) == il.size()){
                        fprintf(file,"%d 1 2 1002 0 %d %d\n", number_of_elements, nodes, last);
                    }else{
                        fprintf(file,"%d 1 2 1002 0 %d %d\n", number_of_elements, nodes, nodes+1);
                    }
            }
        }

        fprintf(file,"$EndElements\n");
        fclose(file);
    //}
}


InspectElementsAlgorithm22::InspectElementsAlgorithm22(Mesh* input_mesh)
: mesh(input_mesh)
{}


void InspectElementsAlgorithm22::compute_intersections(const std::vector< std::vector<ILpair>>& intersection_map_)
{
    DBGMSG("Intersections 2d-2d\n");
    
    FOR_ELEMENTS(mesh, ele) {
    if (ele->dim() == 3)
    {
        unsigned int ele_idx = ele->index();
        // if there are not at least 2 2D elements intersecting 3D element; continue
        if(intersection_map_[ele_idx].size() < 2) continue;
        
        const std::vector<ILpair> &local_map = intersection_map_[ele_idx];
        
        DBGMSG("more than 2 intersections in tetrahedron found\n");
        for(unsigned int i=0; i < local_map.size(); i++)
        {
            //TODO: 1] compute all plucker coords at once
            //TODO: 2] pass plucker coords from 2d-3d
            
            ElementFullIter eleA = mesh->element(local_map[i].first);
            if(eleA->dim() !=2 ) continue;  //skip other dimension intersection
            
            for(unsigned int j=i+1; j < local_map.size(); j++)
            {
                ElementFullIter eleB = mesh->element(local_map[j].first);
                if(eleB->dim() !=2 ) continue;  //skip other dimension intersection
                
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
    ASSERT(eleA->dim() == 2, "Wrong element dimension.");
    ASSERT(eleB->dim() == 2, "Wrong element dimension.");
    ASSERT(eleA->index() != eleB->index(), "Cannot compute intersection of the same elements.");
    
    update_simplex(eleA, triaA_);
    update_simplex(eleB, triaB_);
    
    IntersectionAux<2,2> is(eleA->index(), eleB->index());
    std::vector<unsigned int> prolongation_table;
    
    ComputeIntersection< Simplex<2>, Simplex<2>> CI(triaA_, triaB_);
    CI.init();
    unsigned int n_local_intersection = CI.compute(is, prolongation_table);
    
    if(n_local_intersection > 0)
        intersectionaux_storage22_.push_back(is);
}

void InspectElementsAlgorithm22::update_simplex(const ElementFullIter& element, Simplex< 2 >& simplex)
{
    arma::vec3 *field_of_points[3];
    for(unsigned int i=0; i < 3; i++)
        field_of_points[i]= &(element->node[i]->point());
    simplex.set_simplices(field_of_points);
}


// Declaration of specializations implemented in cpp:
template class InspectElementsAlgorithm<1>;
template class InspectElementsAlgorithm<2>;

} // END namespace
