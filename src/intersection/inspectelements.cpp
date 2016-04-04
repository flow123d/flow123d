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
    intersection_list_.assign(mesh->n_elements(),std::vector<IntersectionAux<dim,3>>());
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
    BIHTree bt(mesh, 20);

    {START_TIMER("Element iteration");
    
    FOR_ELEMENTS(mesh, elm) {
        unsigned int component_ele_idx = elm->index();
        
        if (elm->dim() == dim &&                                // is component element
            !closed_elements[component_ele_idx] &&                    // is not closed yet
            elements_bb[component_ele_idx].intersect(mesh_3D_bb))    // its bounding box intersects 3D mesh bounding box
        {    
            std::vector<unsigned int> searchedElements;
            bt.find_bounding_box(elements_bb[component_ele_idx], searchedElements);

            {START_TIMER("Bounding box element iteration");
            
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
                    
                    DBGMSG("elements %d %d\n",component_ele_idx, bulk_ele_idx);
                    update_simplex(elm, component_simplex); // update component simplex
                    update_simplex(ele_3D, tetrahedron); // update tetrahedron
                    std::vector<unsigned int> prolongation_table;
                    bool found = compute_initial_CI(component_ele_idx, bulk_ele_idx, prolongation_table);

                    // keep the index of the current component element that is being investigated
                    unsigned int current_component_element_idx = component_ele_idx;
                    
                    if(found){
                        
                        last_slave_for_3D_elements[bulk_ele_idx] = component_ele_idx;
                        prolongation_decide(elm, ele_3D, prolongation_table);
                        
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
                        
                        // if component element is closed, do not check other bounding boxes
                        if(closed_elements[component_ele_idx])
                            break;
                    }

                }
            }
            END_TIMER("Bounding box element iteration");}
        }
    }

    END_TIMER("Element iteration");}
    
    // DBG write which elements are closed
    FOR_ELEMENTS(mesh, ele) {
        DBGMSG("Element[%3d] closed: %d\n",ele.index(),(closed_elements[ele.index()] ? 1 : 0));
    }
}
  

 



template<>
void InspectElementsAlgorithm<1>::prolongation_decide(const ElementFullIter& elm, 
                                                      const ElementFullIter& ele_3D,
                                                      std::vector<unsigned int> &prolongation_table)
{
    IntersectionAux<1,3> il = intersection_list_[elm->index()].back();
    // number of IPs that are at vertices of component element
    unsigned int n_ip_vertices = 0;
    
    for(const IntersectionPointAux<1,3> &IP : il.points()) {
        if(IP.dim_A() == 0) {
            n_ip_vertices++;
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
                        DBGMSG("1d prolong %d in %d\n", sousedni_element, ele_3D->index());
                        
                        // Vytvoření průniku bez potřeby počítání
                        IntersectionAux<1,3> il_other(sousedni_element, ele_3D->index());
                        intersection_list_[sousedni_element].push_back(il_other);
                        
                        Prolongation pr = {sousedni_element, ele_3D->index(), (unsigned int)intersection_list_[sousedni_element].size() - 1};
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
                        DBGMSG("3d prolong (node)\n");
                        break;
                
                // IP is on a line of tetrahedron; possible edges are from all connected sides (2)
                case 1: for(unsigned int j=0; j < RefElement<3>::n_sides_per_line; j++)
                            edges.push_back(&(mesh->edges[ele_3D->edge_idx_[RefElement<3>::interact<2,1>(IP.idx_B())[j]]]));
                        DBGMSG("3d prolong (edge)\n");
                        break;
                        
                // IP is on a side of tetrahedron; only possible edge is from the given side
                case 2: edges.push_back(&(mesh->edges[ele_3D->edge_idx_[IP.idx_B()]]));
                        // if edge has only one side, it means it is on the boundary and we cannot prolongate
                        if(edges.back()->n_sides == 1)
                        {
                            Prolongation pr = {elm->index(), undefined_elm_idx_, undefined_elm_idx_};
                            bulk_queue_.push(pr);
                            continue;
                        }
                        DBGMSG("3d prolong (side)\n");
                        break;
                default: ASSERT_LESS(IP.dim_B(),3);
            }
            
            unsigned int n_prolongations = 0;
            for(Edge* edg : edges)
            for(int j=0; j < edg->n_sides;j++) {
                if (edg->side(j)->element() != ele_3D) {
                    unsigned int sousedni_element = edg->side(j)->element()->index();

                    if(last_slave_for_3D_elements[sousedni_element] == undefined_elm_idx_ || 
                        (last_slave_for_3D_elements[sousedni_element] != elm->index() && !intersection_exists(elm->index(),sousedni_element)))
                    {
                        last_slave_for_3D_elements[sousedni_element] = elm->index();
                     
                        DBGMSG("3d prolong %d in %d\n",elm->index(),sousedni_element);
                        
                        // Vytvoření průniku bez potřeby počítání
                        IntersectionAux<1,3> il_other(elm->index(), sousedni_element);
                        intersection_list_[elm->index()].push_back(il_other);
                        
                        Prolongation pr = {elm->index(), sousedni_element, (unsigned int)intersection_list_[elm->index()].size() - 1};
                        bulk_queue_.push(pr);
                        n_prolongations++;
                    }
                }   
            }
            
            DBGMSG("cover: %d %d\n", il.size(), n_prolongations);
            // if there are no sides of any edge that we can continue to prolongate over,
            // it means we are at the boundary and cannot prolongate further
            if(IP.is_pathologic() && n_prolongations == 0)
            {
                Prolongation pr = {elm->index(), undefined_elm_idx_, undefined_elm_idx_};
                bulk_queue_.push(pr);
            }
        }  
    }
    
    if(n_ip_vertices == il.size()) closed_elements[elm->index()] = true;
}


template<>
void InspectElementsAlgorithm<2>::prolongation_decide(const ElementFullIter& elm, 
                                                      const ElementFullIter& ele_3D,
                                                      std::vector<unsigned int> &prolongation_table
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

            SideIter elm_side = elm->side(side);
            Edge *edg = elm_side->edge();

            for(int j=0; j < edg->n_sides;j++) {
                SideIter other_side=edg->side(j);
                if (other_side != elm_side) {
                    unsigned int sousedni_element = other_side->element()->index(); // 2D element

                    DBGMSG("2d sousedni_element %d\n", sousedni_element);

                        if(!intersection_exists(sousedni_element,ele_3D->index())){

                            DBGMSG("2d prolong\n");
                            // Vytvoření průniku bez potřeby počítání
                            IntersectionAux<2,3> il_other(sousedni_element, ele_3D->index());
                            intersection_list_[sousedni_element].push_back(il_other);

                            Prolongation pr = {sousedni_element, ele_3D->index(), (unsigned int)intersection_list_[sousedni_element].size() - 1};
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

                    // TODO:
                    // - rename it, describe it and test that it is really useful !!
                    // how it probably works: 
                    // - if index undefined then no intersection has been computed for the 3D element
                    // - if not undefined check the index of actual 2D element 
                    //   (most probable case, that we are looking along 2D element back to 3D element, which we have just computed)
                    // - if it is another 2D element, then go through all found intersections of the 3D element and test it..
                    
                    if(last_slave_for_3D_elements[sousedni_element] == undefined_elm_idx_ || 
                        (last_slave_for_3D_elements[sousedni_element] != elm->index() && !intersection_exists(elm->index(),sousedni_element))){
                        
                        last_slave_for_3D_elements[sousedni_element] = elm->index();

                        DBGMSG("3d prolong\n");
                        
                        // Vytvoření průniku bez potřeby počítání
                        IntersectionAux<2,3> il_other(elm->index(), sousedni_element);
                        intersection_list_[elm->index()].push_back(il_other);

                        Prolongation pr = {elm->index(), sousedni_element, (unsigned int)intersection_list_[elm->index()].size() - 1};
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
    
    ComputeIntersection<Simplex<dim>, Simplex<3>> CI(component_simplex, tetrahedron);
    CI.init();
    CI.compute(is, prolongation_table);

    if(is.size() > dim){
//         for(unsigned int j=1; j < is.size(); j++) 
//             cout << is[j];
        
        prolongation_decide(elm, ele_3D, prolongation_table);
        n_intersections_++;
    }
}


 
 
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
 

template<unsigned int dim>
void InspectElements::compute_intersections(std::vector<IntersectionLocal<dim,3>> &storage)
{
    InspectElementsAlgorithm<dim> iea(mesh);  
    iea.compute_intersections();
  
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
    
}
 
void InspectElements::compute_intersections()
{
    intersection_map_.resize(mesh->n_elements());
    
    compute_intersections<1>(intersection_storage13_);
    compute_intersections<2>(intersection_storage23_);
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
            }
        }

        unsigned int number_of_elements = mesh->n_elements();
        unsigned int nodes = mesh->n_nodes();

        for(unsigned int j = 0; j < intersection_storage13_.size();j++){
            IntersectionLocal<1,3> il = intersection_storage13_[j];
            number_of_elements++;
            nodes++;
            if(il.size() == 1){
                fprintf(file,"%d 1 2 18 7 %d %d\n", number_of_elements, nodes, nodes);
            }else if(il.size() == 2){
                fprintf(file,"%d 1 2 18 7 %d %d\n", number_of_elements, nodes, nodes+1);
                nodes++;
            }
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
