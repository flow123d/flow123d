/*
 * inspectelements.cpp
 *
 *  Created on: 13.4.2014
 *      Author: viktor, pe, jb
 */

#include "inspectelements.h"
#include "inspect_elements_algorithm.h"
#include "intersectionpoint.h"
#include "intersectionaux.h"
#include "intersection_local.h"

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"

namespace computeintersection {

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
//     iea.compute_intersections_BIHtree();
//     iea.compute_intersections_BB();
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
//                     IntersectionLocal<dim,3>* il = 
//                         static_cast<IntersectionLocal<dim,3>*> (intersection_map_[idx][j].second);
//                     cout << il;
//                     for(IntersectionPoint<dim,3> &ip : il->points())
//                         ip.coords(mesh->element(idx)).print(cout);

                    // create map for bulk element
                    intersection_map_[bulk_idx].push_back(
                                                std::make_pair(
                                                    idx, 
                                                    &(storage.back())
                                                ));
                }
        }
    }
    END_TIMER("Intersection into storage");
    
    for(IntersectionLocal<2,3> &is : intersection_storage23_) {
        xprintf(Msg, "comp-bulk: %4d %4d\n", is.component_ele_idx(), is.bulk_ele_idx());
    }
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
    
    if(d & (IntersectionType::d23 | IntersectionType::d22)){
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

        unsigned int idx = 1;
        FOR_NODES(mesh, nod){
            arma::vec3 _nod = nod->point();
            fprintf(file,"%d %.16f %.16f %.16f\n", idx, _nod[0], _nod[1], _nod[2]);
            idx++;
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

                fprintf(file,"%d %.16f %.16f %.16f\n", number_of_nodes, global[0], global[1], global[2]);
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

        unsigned int idx = 1;
        FOR_NODES(mesh, nod){
            arma::vec3 _nod = nod->point();
            fprintf(file,"%d %.16f %.16f %.16f\n", idx, _nod[0], _nod[1], _nod[2]);
            idx++;
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
                    fprintf(file,"%d %.16f %.16f %.16f\n", number_of_nodes, global[0], global[1], global[2]);
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

                    if((k+1) == il.size() && il.size()){
                        fprintf(file,"%d 1 2 1002 0 %d %d\n", number_of_elements, nodes, last);
                    }else{
                        fprintf(file,"%d 1 2 1002 0 %d %d\n", number_of_elements, nodes, nodes+1);
                    }
                    //if(il.size() < 3) break;
            }
        }

        fprintf(file,"$EndElements\n");
        fclose(file);
    //}
}

} // END namespace
