/*
 * inspectelements.cpp
 *
 *  Created on: 13.4.2014
 *      Author: viktor
 */

#include "inspectelements.h"
#include "intersectionpoint.h"
#include "trace_algorithm.h"
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
: component_element_idx_(-1), mesh(_mesh)
{
    intersection_list_.assign(mesh->n_elements(),std::vector<IntersectionAux<dim,3>>());
}

template<unsigned int dim>   
InspectElementsAlgorithm<dim>::~InspectElementsAlgorithm()
{}


    
template<unsigned int dim>
void InspectElementsAlgorithm<dim>::init()
{
    START_TIMER("Intersection initialization");
    last_slave_for_3D_elements.assign(mesh->n_elements(), -1);
    closed_elements.assign(mesh->n_elements(), false);

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


template<>
void InspectElementsAlgorithm<1>::trace(IntersectionAux<1,3> &intersection)
{}
 
template<>
void InspectElementsAlgorithm<2>::trace(IntersectionAux<2,3> &intersection)
{
    prolongation_table_.clear();
    Tracing::trace_polygon(prolongation_table_, intersection);
} 
 
template<unsigned int dim> 
bool InspectElementsAlgorithm<dim>::compute_initial_CI(const ElementFullIter& elm, const ElementFullIter& ele_3D)
{
    IntersectionAux<dim,3> is(elm->index(), ele_3D->index());
    START_TIMER("Compute intersection");
    ComputeIntersection<Simplex<dim>, Simplex<3>> CI(component_simplex, tetrahedron);
    CI.init();
    CI.compute(is.points());
    END_TIMER("Compute intersection");
    
    if(is.points().size() > 0) {
        
        trace(is);
        intersection_list_[elm->index()].push_back(is);
        return true;
    }
    else return false;
}

template<unsigned int dim>
bool InspectElementsAlgorithm<dim>::intersection_exists(unsigned int component_element_idx, unsigned int elm_3D_idx) 
{
    bool found = false;

    for(unsigned int i = 0; i < intersection_list_[component_element_idx].size();i++){

        if(intersection_list_[component_element_idx][i].bulk_ele_idx() == elm_3D_idx){
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
    
    BIHTree bt(mesh, 20);

    {START_TIMER("Element iteration");
    
    FOR_ELEMENTS(mesh, elm) {
        if (elm->dim() == dim &&                                // is 2D element
            !closed_elements[elm.index()] &&                    // is not closed yet
            elements_bb[elm->index()].intersect(mesh_3D_bb))    // its bounding box intersects 3D mesh bounding box
        {    
            update_simplex(elm, component_simplex); // update component simplex
            std::vector<unsigned int> searchedElements;
            bt.find_bounding_box(elements_bb[elm->index()], searchedElements);

            {START_TIMER("Bounding box element iteration");
            
            for (std::vector<unsigned int>::iterator it = searchedElements.begin(); it!=searchedElements.end(); it++)
            {
                int idx = *it;
                ElementFullIter ele_3D = mesh->element(idx);

                if (ele_3D->dim() == 3 && last_slave_for_3D_elements[ele_3D->index()] != (int)(elm->index())) {
                    
                    DBGMSG("elements %d %d\n",elm->index(), ele_3D->index());
                    update_simplex(ele_3D, tetrahedron); // update tetrahedron
                    bool found = compute_initial_CI(elm, ele_3D);

                    //TODO: how does prolongation work when there are 1 or 2 IPs?
                    if(found){
                        
                        last_slave_for_3D_elements[ele_3D->index()] = elm.index();

                        // prolongation:
                        
                        prolongation_decide(elm, ele_3D);

                        // - find first intersection polygon
                        // - fill both prolongation queue in 2D and 3D
                        // - clear prolongation queue:
                        //      - clear 3D prolongation queue:
                        //          - compute CI
                        //          - fill both prolongation queue 2D and 3D
                        //          - repeat until 3D queue is empty
                        //      - we can close the 2D element which we started with
                        //      - clear 2D prolongation queue:
                        //          - compute CI
                        //          - fill both prolongation queue 2D and 3D
                        //          - repeat until 2D queue is empty
                        //      - emptying 2D queue could fill 3D queue again, so repeat
                        while(1){
                            while(!bulk_queue_.empty()){
                                Prolongation pr = bulk_queue_.front();
                                DBGMSG("Prolongation queue compute in 3d ele_idx %d.\n",pr.elm_3D_idx);
                                prolongate(pr);
                                bulk_queue_.pop();

                            }

                            if(component_element_idx_ >= 0) closed_elements[component_element_idx_] = true;

                            component_element_idx_ = -1;
                            
                            if(!component_queue_.empty()){
                                Prolongation pr = component_queue_.front();
                                DBGMSG("Prolongation queue compute in 2d ele_idx %d.\n",pr.component_elm_idx);
                                prolongate(pr);
                                component_queue_.pop();

                            }

                                
                            if(component_queue_.empty() && bulk_queue_.empty()){
                                if(component_element_idx_ >= 0) closed_elements[component_element_idx_] = true;
                                    
                                component_element_idx_ = -1;
                                break;
                            }

                        }
                        break; // ukončí procházení dalších bounding boxů
                    }

                }
            }
            // Prošlo se celé pole sousedním bounding boxů, pokud nevznikl průnik, může se trojúhelník nastavit jako uzavřený
            closed_elements[elm.index()] = true;
            END_TIMER("Bounding box element iteration");}
        }
    }

    END_TIMER("Element iteration");}
}
  

 



template<>
void InspectElementsAlgorithm<1>::prolongation_decide(const ElementFullIter& elm, const ElementFullIter& ele_3D)
{
    IntersectionAux<1,3> il = intersection_list_[elm.index()].back();
    for(const IntersectionPoint<1,3> &IP : il.points()) {
        if(IP.dim_A() == 0) { 
            DBGMSG("1D end\n");
            // if IP is end of the 1D element
            // prolongate 1D element as long as it creates prolongation point on the side of tetrahedron
            SideIter elm_side = elm->side(IP.idx_A());  // side of 1D element is vertex
            Edge *edg = elm_side->edge();
            for(int j=0; j < edg->n_sides;j++) {
                
                SideIter other_side=edg->side(j);
                if (other_side != elm_side) {

                    unsigned int sousedni_element = other_side->element()->index();
                    if(!intersection_exists(sousedni_element,ele_3D->index())){
                    //if(!closed_elements[sousedni_element]){
                        DBGMSG("1d prolong %d in %d\n", sousedni_element, ele_3D->index());
                        
                        // Vytvoření průniku bez potřeby počítání
                        IntersectionAux<1,3> il_other(sousedni_element, ele_3D->index());
                        intersection_list_[sousedni_element].push_back(il_other);
                        
                        Prolongation pr = {sousedni_element, ele_3D->index(), intersection_list_[sousedni_element].size() - 1};
                        component_queue_.push(pr);
                    }
                }
            }
        }else{
            std::vector<Edge*> edges;

            switch (IP.dim_B())
            {
                // IP is at a node of tetrahedron; possible edges are from all connected sides (3)
                case 0: for(unsigned int j=0; j < RefElement<3>::n_sides_per_node; j++)
                            edges.push_back(&(mesh->edges[ele_3D->edge_idx_[RefElement<3>::interact<2,0>(IP.idx_B())[j]]]));
                        break;
                
                // IP is on a line of tetrahedron; possible edges are from all connected sides (2)
                case 1: for(unsigned int j=0; j < RefElement<3>::n_sides_per_line; j++)
                            edges.push_back(&(mesh->edges[ele_3D->edge_idx_[RefElement<3>::interact<2,1>(IP.idx_B())[j]]]));
                        break;
                        
                // IP is on a side of tetrahedron; only possible edge is from the given side
                case 2: edges.push_back(&(mesh->edges[ele_3D->edge_idx_[IP.idx_B()]]));
                        break;
                default: ASSERT_LESS(IP.dim_B(),3);
            }
            
            for(Edge* edg : edges)
            for(int j=0; j < edg->n_sides;j++) {
                if (edg->side(j)->element() != ele_3D) {
                    unsigned int sousedni_element = edg->side(j)->element()->index();

                    if(last_slave_for_3D_elements[sousedni_element] == -1 || 
                        (last_slave_for_3D_elements[sousedni_element] != (int)elm->index() && !intersection_exists(elm->index(),sousedni_element)))
                    {
                        last_slave_for_3D_elements[sousedni_element] = elm->index();
                     
                        DBGMSG("3d prolong\n");
                        
                        // Vytvoření průniku bez potřeby počítání
                        IntersectionAux<1,3> il_other(elm.index(), sousedni_element);
                        intersection_list_[elm.index()].push_back(il_other);
                        
                        Prolongation pr = {elm->index(), sousedni_element, intersection_list_[elm.index()].size() - 1};
                        bulk_queue_.push(pr);
                    }
                }   
            }
        }  
    }
}


template<>
void InspectElementsAlgorithm<2>::prolongation_decide(const ElementFullIter& elm, const ElementFullIter& ele_3D)
{
    for(unsigned int i = 0; i < prolongation_table_.size();i++){

        unsigned int side;
        bool is_triangle_side = true;

        if(prolongation_table_[i] >= 4){
            side = prolongation_table_[i] - 4;
        }else{
            side = prolongation_table_[i];
            is_triangle_side = false;
        }

        DBGMSG("prolongation table: %d %d\n", side, is_triangle_side);

        if(is_triangle_side){
            // prolongation through the triangle side

            SideIter elm_side = elm->side(side);


            Edge *edg = elm_side->edge();

            for(int j=0; j < edg->n_sides;j++) {
                SideIter other_side=edg->side(j);
                if (other_side != elm_side) {
                    unsigned int sousedni_element = other_side->element()->index(); // 2D element

                    DBGMSG("2d sousedni_element %d\n", sousedni_element);
                            //xprintf(Msg, "Naleznut sousedni element elementu(3030) s ctyrstenem(%d)\n", ele->index());
                            //xprintf(Msg, "\t\t Idx původního elementu a jeho hrany(%d,%d) - Idx sousedního elementu a jeho hrany(%d,%d)\n",elm->index(),side,other_side->element()->index(),other_side->el_idx());

                        if(!intersection_exists(sousedni_element,ele_3D->index())){
                            //flag_for_3D_elements[ele->index()] = sousedni_element;
                            DBGMSG("2d prolong\n");
                            // Vytvoření průniku bez potřeby počítání
                            IntersectionAux<2,3> il_other(sousedni_element, ele_3D->index());
                            intersection_list_[sousedni_element].push_back(il_other);

                            Prolongation pr = {sousedni_element, ele_3D->index(), intersection_list_[sousedni_element].size() - 1};
                            component_queue_.push(pr);
                        }
                }
            }

        }else{
            // prolongation through the tetrahedron side
            
            SideIter elm_side = ele_3D->side(side);
            Edge *edg = elm_side->edge();

            for(int j=0; j < edg->n_sides;j++) {
                SideIter other_side=edg->side(j);
                if (other_side != elm_side) {
                    //

                    unsigned int sousedni_element = other_side->element()->index();

                    DBGMSG("3d sousedni_element %d\n", sousedni_element);

                    // TODO: flag_for_3D_elements seems to be optimalisation:
                    // - rename it, describe it and test that it is really useful !!
                    // how it probably works: 
                    // - if "-1" then no intersection has been computed for the 3D element
                    // - if not "-1" check the index of actual 2D element 
                    //   (most probable case, that we are looking along 2D element back to 3D element, which we have just computed)
                    // - if it is another 2D element, then go through all found intersections of the 3D element and test it..
                    
                    if(last_slave_for_3D_elements[sousedni_element] == -1 || 
                        (last_slave_for_3D_elements[sousedni_element] != (int)elm->index() && !intersection_exists(elm->index(),sousedni_element))){
                        
                        last_slave_for_3D_elements[sousedni_element] = elm->index();
                        // Jedná se o vnitřní čtyřstěny v trojúhelníku

                        DBGMSG("3d prolong\n");
                        
                        // Vytvoření průniku bez potřeby počítání
                        IntersectionAux<2,3> il_other(elm.index(), sousedni_element);
                        intersection_list_[elm.index()].push_back(il_other);

                        Prolongation pr = {elm.index(), sousedni_element, intersection_list_[elm.index()].size() - 1};
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
    
    update_simplex(elm, component_simplex);
    update_simplex(ele_3D, tetrahedron);
    
    component_element_idx_ = pr.component_elm_idx;

    IntersectionAux<dim,3> &is = intersection_list_[pr.component_elm_idx][pr.dictionary_idx];
    
    ComputeIntersection<Simplex<dim>, Simplex<3>> CI(component_simplex, tetrahedron);
    CI.init();
    CI.compute(is.points());

    if(is.size() > dim){
//         for(unsigned int j=1; j < is.size(); j++) 
//             cout << is[j];
        
        trace(is);
        prolongation_decide(elm, ele_3D);
    }
}


// template<> double InspectElementsAlgorithm<1>::measure() 
// {
//     double subtotal = 0.0;
// 
//     for(unsigned int i = 0; i < intersection_list_.size(); i++){
//         DBGMSG("intersection_list_[i].size() = %d\n",intersection_list_[i].size());
//         for(unsigned int j = 0; j < intersection_list_[i].size();j++){
//             //if(intersection_list_[i][j].size() < 2) continue;
//             
//             ElementFullIter ele = mesh->element(intersection_list_[i][j].component_ele_idx());            
//             double t1d_length = ele->measure();
//             double local_length = intersection_list_[i][j].compute_measure();
//             
//             if(intersection_list_[i][j].size() == 2)
//             {
//             arma::vec3 from = intersection_list_[i][j][0].coords(ele);
//             arma::vec3 to = intersection_list_[i][j][1].coords(ele);
//             DBGMSG("sublength from [%f %f %f] to [%f %f %f] = %f\n",
//                    from[0], from[1], from[2], 
//                    to[0], to[1], to[2],
//                    local_length*t1d_length);
//             }
//             subtotal += local_length*t1d_length;
//         }
//     }
//     return subtotal;
// }
// 
// template<> double InspectElementsAlgorithm<2>::measure() 
// {
//     double subtotal = 0.0;
// 
//     for(unsigned int i = 0; i < intersection_list_.size(); i++){
//         for(unsigned int j = 0; j < intersection_list_[i].size();j++){
//             double t2dArea = mesh->element(intersection_list_[i][j].component_ele_idx())->measure();
//             double localArea = intersection_list_[i][j].compute_measure();
//             subtotal += 2*localArea*t2dArea;
//         }
//     }
//     return subtotal;
// } 
 
 
 
 
 
 
 
InspectElements::InspectElements(Mesh* mesh)
: mesh(mesh)
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
 
void InspectElements::compute_intersections()
{
    InspectElementsAlgorithm<1> iea_13(mesh);
    InspectElementsAlgorithm<2> iea_23(mesh);
    
    iea_13.compute_intersections();
    iea_23.compute_intersections();
    
    intersection_map_.resize(mesh->n_elements());
    intersection_storage13_.reserve(100);
    intersection_storage23_.reserve(100);
    
    FOR_ELEMENTS(mesh, elm) {
        unsigned int idx = elm->index(); 
        unsigned int bulk_idx;
        
        switch (elm->dim()) {   
            case 1: {
                intersection_map_[idx].resize(iea_13.intersection_list_[idx].size());
                for(unsigned int j = 0; j < iea_13.intersection_list_[idx].size(); j++){
                    bulk_idx = iea_13.intersection_list_[idx][j].bulk_ele_idx();
//                     DBGMSG("cidx %d bidx %d: n=%d\n",idx, bulk_idx, iea_13.intersection_list_[idx][j].size());
                    intersection_storage13_.push_back(IntersectionLocal<1,3>(iea_13.intersection_list_[idx][j]));
                    
                    // create map for component element
                    intersection_map_[idx][j] = std::make_pair(
                                                    bulk_idx, 
                                                    &(intersection_storage13_.back())
                                                );
//                  // write down intersections
//                     IntersectionLocal<1,3>* il13 = 
//                         static_cast<IntersectionLocal<1,3>*> (intersection_map_[idx][j].second);
//                     cout << &(intersection_storage13_.back()) << "  " << il13 << *il13;
//                     for(IntersectionPointX<1,3> &ip : il13->points())
//                         ip.coords(mesh->element(idx)).print(cout);

                    // create map for bulk element
                    intersection_map_[bulk_idx].push_back(
                                                std::make_pair(
                                                    idx, 
                                                    &(intersection_storage13_.back())
                                                ));
                }
                break;
            }
            case 2: {
                intersection_map_[idx].resize(iea_23.intersection_list_[idx].size());
                for(unsigned int j = 0; j < iea_23.intersection_list_[idx].size(); j++){
                    bulk_idx = iea_23.intersection_list_[idx][j].bulk_ele_idx();
                    intersection_storage23_.push_back(IntersectionLocal<2,3>(iea_23.intersection_list_[idx][j]));
                    
                    // create map for component element
                    intersection_map_[idx][j] = std::make_pair(
                                                    bulk_idx, 
                                                    (IntersectionLocalBase*) &(intersection_storage23_.back())
                                                );
                    // create map for bulk element
                    intersection_map_[bulk_idx].push_back(
                                                std::make_pair(
                                                    idx, 
                                                    (IntersectionLocalBase*) &(intersection_storage23_.back())
                                                ));
                }
                break;
            }
            default: break;
            
        }
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
            ElementFullIter el3D = mesh->element(il.bulk_ele_idx());
            
            for(unsigned int k = 0; k < il.size();k++){
                number_of_nodes++;
                IntersectionPointX<1,3> IP13 = il[k];
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
        fprintf(file,"%d\n", (intersection_storage13_.size() + mesh->n_elements()) );

        FOR_ELEMENTS(mesh, elee){
            // 1 4 2 30 26 1 2 3 4
            // 2 2 2 2 36 5 6 7
            if(elee->dim() == 3){
                int id1 = mesh->node_vector.index(elee->node[0]) + 1;
                int id2 = mesh->node_vector.index(elee->node[1]) + 1;
                int id3 = mesh->node_vector.index(elee->node[2]) + 1;
                int id4 = mesh->node_vector.index(elee->node[3]) + 1;

                fprintf(file,"%d 4 2 30 26 %d %d %d %d\n", elee.id(), id1, id2, id3, id4);
            }else if(elee->dim() == 2){
                int id1 = mesh->node_vector.index(elee->node[0]) + 1;
                int id2 = mesh->node_vector.index(elee->node[1]) + 1;
                int id3 = mesh->node_vector.index(elee->node[2]) + 1;
                fprintf(file,"%d 2 2 2 36 %d %d %d\n", elee.id(), id1, id2, id3);

            }else if(elee->dim() == 1){
                int id1 = mesh->node_vector.index(elee->node[0]) + 1;
                int id2 = mesh->node_vector.index(elee->node[1]) + 1;
                fprintf(file,"%d 1 2 14 16 %d %d\n",elee.id(), id1, id2);
                //xprintf(Msg, "Missing implementation of printing 1D element to a file\n");
            }
        }

        unsigned int number_of_elements = mesh->n_elements();
        //unsigned int last = 0;
        unsigned int nodes = mesh->n_nodes();

        for(unsigned int j = 0; j < intersection_storage13_.size();j++){
            IntersectionLocal<1,3> il = intersection_storage13_[j];
            //for(unsigned int k = 0; k < il.size();k++){

            number_of_elements++;
            nodes++;
            if(il.size() == 1){
                fprintf(file,"%d 1 2 18 7 %d %d\n", number_of_elements, nodes, nodes);
            }else if(il.size() == 2){
                fprintf(file,"%d 1 2 18 7 %d %d\n", number_of_elements, nodes, nodes+1);
                nodes++;
            }
            //}
        }

        fprintf(file,"$EndElements\n");
        fclose(file);
}

void InspectElements::print_mesh_to_file_23(string name)
{
    for(unsigned int i = 0; i < 2;i++){
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
            ElementFullIter el3D = mesh->element(il.bulk_ele_idx());
            
            for(unsigned int k = 0; k < intersection_storage23_[j].size();k++){

                    number_of_nodes++;
                    IntersectionPointX<2,3> IP23 = il[k];
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
            // 1 4 2 30 26 1 2 3 4
            // 2 2 2 2 36 5 6 7
            if(elee->dim() == 3){
                int id1 = mesh->node_vector.index(elee->node[0]) + 1;
                int id2 = mesh->node_vector.index(elee->node[1]) + 1;
                int id3 = mesh->node_vector.index(elee->node[2]) + 1;
                int id4 = mesh->node_vector.index(elee->node[3]) + 1;

                fprintf(file,"%d 4 2 30 26 %d %d %d %d\n", elee.id(), id1, id2, id3, id4);
            }else if(elee->dim() == 2){
                int id1 = mesh->node_vector.index(elee->node[0]) + 1;
                int id2 = mesh->node_vector.index(elee->node[1]) + 1;
                int id3 = mesh->node_vector.index(elee->node[2]) + 1;
                fprintf(file,"%d 2 2 2 36 %d %d %d\n", elee.id(), id1, id2, id3);

            }else{
                int id1 = mesh->node_vector.index(elee->node[0]) + 1;
                int id2 = mesh->node_vector.index(elee->node[1]) + 1;
                fprintf(file,"%d 1 2 14 16 %d %d\n",elee.id(), id1, id2);
                //xprintf(Msg, "Missing implementation of printing 1D element to a file\n");
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
                        fprintf(file,"%d 1 2 18 7 %d %d\n", number_of_elements, nodes, last);
                    }else{
                        fprintf(file,"%d 1 2 18 7 %d %d\n", number_of_elements, nodes, nodes+1);
                    }
            }
        }

        fprintf(file,"$EndElements\n");
        fclose(file);
    }
}



// Declaration of specializations implemented in cpp:
template class InspectElementsAlgorithm<1>;
template class InspectElementsAlgorithm<2>;

} // END namespace
