/*
 * MixedMeshIntersections.cpp
 *
 *  Created on: 13.4.2014
 *      Author: viktor, pe, jb
 */

#include "inspect_elements_algorithm.hh"
#include "intersection_point_aux.hh"
#include "intersection_aux.hh"
#include "intersection_local.hh"

#include "system/global_defs.h"
#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "mesh/ref_element.hh"
#include "mixed_mesh_intersections.hh"
#include "mesh/bih_tree.hh"
#include "mesh/accessors.hh"
#include "mesh/node_accessor.hh"
#include "mesh/range_wrapper.hh"


MixedMeshIntersections::MixedMeshIntersections(Mesh* mesh)
: mesh(mesh), algorithm13_(mesh), algorithm23_(mesh), algorithm22_(mesh), algorithm12_(mesh)
{}

MixedMeshIntersections::~MixedMeshIntersections()
{}


unsigned int MixedMeshIntersections::number_of_components(unsigned int dim)
{
    ASSERT(dim < 3);
//     if(dim == 1) return algorithm13_.component_counter_;
    if(dim == 2) return algorithm22_.component_counter_;
    else xprintf(Err, "Not implemented for dim %d\n.",dim);
    return 0;
}

 
double MixedMeshIntersections::measure_13()
{
    double subtotal = 0.0;

    for(unsigned int i = 0; i < intersection_storage13_.size(); i++){
    	ElementAccessor<3> ele = mesh->element_accessor( intersection_storage13_[i].component_ele_idx() );
        double t1d_length = ele.measure();
        double local_length = intersection_storage13_[i].compute_measure();
        
        if(intersection_storage13_[i].size() == 2)
        {
        arma::vec3 from = intersection_storage13_[i][0].coords(ele);
        arma::vec3 to = intersection_storage13_[i][1].coords(ele);
//         DebugOut().fmt("sublength from [{} {} {}] to [{} {} {}] = %f\n",
//                from[0], from[1], from[2],
//                to[0], to[1], to[2],
//                local_length*t1d_length);
        }
        subtotal += local_length*t1d_length;
    }
    return subtotal;
}


double MixedMeshIntersections::measure_23()
{
    double subtotal = 0.0;

    for(unsigned int i = 0; i < intersection_storage23_.size(); i++){
            double t2dArea = mesh->element_accessor( intersection_storage23_[i].component_ele_idx() ).measure();
            double localArea = intersection_storage23_[i].compute_measure();
            subtotal += 2*localArea*t2dArea;
        }
    return subtotal;
} 


double MixedMeshIntersections::measure_22()
{
    double subtotal = 0.0;
    double val;
    
    for(unsigned int i = 0; i < intersection_storage22_.size(); i++){
        if(intersection_storage22_[i].size() > 1)
        {
        	ElementAccessor<3> eleA = mesh->element_accessor( intersection_storage22_[i].component_ele_idx() );
//             ElementAccessor<3> eleB = mesh->element_accessor( intersection_storage22_[i].bulk_ele_idx() );
            
            arma::vec3 from = intersection_storage22_[i][0].coords(eleA);
            arma::vec3 to = intersection_storage22_[i][1].coords(eleA);
            val = arma::norm(from - to, 2);

//             DebugOut().fmt("{}--{}:sublength from [{} {} {}] to [{} {} {}] = {}\n",
//                eleA.idx(), eleB.idx(),
//                from[0], from[1], from[2],
//                to[0], to[1], to[2],
//                val);
            subtotal += val;
        }
    }
    return subtotal;
}

 
template<uint dim_A, uint dim_B>
void MixedMeshIntersections::store_intersection(std::vector<IntersectionLocal<dim_A, dim_B>> &storage, IntersectionAux<dim_A, dim_B> &isec_aux) {
    //unsigned int ele_a_idx = isec_aux.component_ele_idx();
    //unsigned int ele_b_idx = isec_aux.bulk_ele_idx();

    //WARNING: 
    // - not all algorithms uses this function (e.g. 2d2d pushes directly into storage)
    // - we cannot throw away isec of zero measure generaly (e.g. 1d2d, 1d3d)
    // - is it not better to test number of IPs according to dimensions?
    
    IntersectionLocal<dim_A, dim_B> isec(isec_aux);
    //if ( (! (dim_A==1 && dim_B==2)) && isec.compute_measure() < 1e-14) return;

    storage.push_back(isec);
    /*
    element_intersections_[ele_a_idx].push_back(
            std::make_pair(ele_b_idx, &(storage.back())) );

    element_intersections_[ele_b_idx].push_back(
            std::make_pair(ele_a_idx, &(storage.back())) );
    */

}

template<uint dim_A, uint dim_B>
void MixedMeshIntersections::append_to_index( std::vector<IntersectionLocal<dim_A, dim_B>> &storage)
{
    for(auto &isec : storage) {

        unsigned int ele_a_idx = isec.component_ele_idx();
        unsigned int ele_b_idx = isec.bulk_ele_idx();
        ASSERT_EQ_DBG(mesh->element_accessor(ele_a_idx)->dim(), dim_A)(ele_a_idx);
        ASSERT_EQ_DBG(mesh->element_accessor(ele_b_idx)->dim(), dim_B)(ele_b_idx);
        element_intersections_[ele_a_idx].push_back(
                    std::make_pair(ele_b_idx, &(isec)) );


        if (dim_B==3) {
            // necessary for 2d-2d intersections
            element_intersections_[ele_b_idx].push_back(
                        std::make_pair(ele_a_idx, &(isec)) );

        }
/*
        element_intersections_[ele_b_idx].push_back(
                std::make_pair(ele_a_idx, &(isec)) );*/
    }
}


template<unsigned int dim>
void MixedMeshIntersections::compute_intersections(InspectElementsAlgorithm< dim >& iea,
                                            std::vector< IntersectionLocal<dim,3>>& storage)
{
    START_TIMER("Intersection algorithm");

    Mesh::IntersectionSearch is = mesh->get_intersection_search();
    switch(is){
        case Mesh::BIHsearch: iea.compute_intersections(mesh->get_bih_tree()); break;
        case Mesh::BIHonly:   iea.compute_intersections_BIHtree(mesh->get_bih_tree()); break;
        case Mesh::BBsearch:  iea.compute_intersections_BB(); break;
        default: ASSERT(0).error("Unsupported search algorithm.");
    }
    
    END_TIMER("Intersection algorithm");
    
    START_TIMER("Intersection into storage");
    storage.reserve(iea.n_intersections_);
    
    for (auto elm : mesh->elements_range()) {
        unsigned int idx = elm.idx();
        
        if(elm->dim() == dim)
        {
//                 intersection_map_[idx].resize(iea.intersection_list_[idx].size());
                element_intersections_[idx].reserve(iea.intersection_list_[idx].size());
                for(unsigned int j = 0; j < iea.intersection_list_[idx].size(); j++){
                    
                    // skip zero intersections (are made in iea.prolongate())
                    if(iea.intersection_list_[idx][j].size() == 0) continue;
                    store_intersection(storage, iea.intersection_list_[idx][j]);                }
        }
    }
    END_TIMER("Intersection into storage");
    
//     for(IntersectionLocal<2,3> &is : intersection_storage23_) {
//         DebugOut().fmt("comp-bulk: {} {}\n", is.component_ele_idx(), is.bulk_ele_idx());
//     }
}

void MixedMeshIntersections::compute_intersections_22(vector< IntersectionLocal< 2, 2 > >& storage)
{
    START_TIMER("Intersection algorithm");
    algorithm22_.compute_intersections(element_intersections_, storage);
    END_TIMER("Intersection algorithm");
    
//     START_TIMER("Intersection into storage");
// 
//     storage.reserve(algorithm22_.intersectionaux_storage22_.size());
//     
//     for(IntersectionAux<2,2> &is : algorithm22_.intersectionaux_storage22_) {
//         unsigned int triaA_idx = is.component_ele_idx();
//         unsigned int triaB_idx = is.bulk_ele_idx();
// 
//         //HACK: 'skip flag' move this check into algorithm12_.compute_intersections()
//         bool skip = false;
//         for(unsigned int i=0; i<element_intersections_[triaA_idx].size(); i++)
//         {
//             if(element_intersections_[triaA_idx][i].first == triaB_idx)
//                 skip = true;
//         }
//         if(! skip) {
//             storage.push_back(IntersectionLocal<2,2>(is));
//             element_intersections_[triaA_idx].push_back(std::make_pair(
//                                                         triaB_idx,
//                                                         &(storage.back())
//                                                     ));
//             element_intersections_[triaB_idx].push_back(std::make_pair(
//                                                         triaA_idx,
//                                                         &(storage.back())
//                                                     ));
//         
// //             DebugOut().fmt("2D-2D intersection [{} - {}]:\n",
// //                            mesh->element_accessor(is.component_ele_idx()).idx(),
// //                            mesh->element_accessor(is.bulk_ele_idx()).idx());
// //             for(const IntersectionPointAux<2,2>& ip : is.points()) {
// //                 DebugOut() << ip;
// //                 auto p = ip.coords(mesh->element(is.component_ele_idx()));
// //                 DebugOut() << "[" << p[0] << " " << p[1] << " " << p[2] << "]\n";
// //             }
//         }
//     }
//     DBGVAR(algorithm22_.intersectionaux_storage22_.size());
// 
//     END_TIMER("Intersection into storage");
}

void MixedMeshIntersections::compute_intersections_12_3(vector< IntersectionLocal< 1, 2 > >& storage)
{
    storage.reserve(intersection_storage13_.size());
    algorithm12_.compute_intersections_3(element_intersections_, storage);
    storage.shrink_to_fit();
    
//     START_TIMER("Intersection into storage");
//     storage.reserve(algorithm12_.intersectionaux_storage12_.size());
    
//     for(IntersectionAux<1,2> &is : algorithm12_.intersectionaux_storage12_) {
//         unsigned int abscissa_idx = is.component_ele_idx();
//         unsigned int triangle_idx = is.bulk_ele_idx();
// 
//         //HACK: 'skip flag' move this check into algorithm12_.compute_intersections()
//         bool skip = false;
//         for(unsigned int i=0; i<intersection_map_[abscissa_idx].size(); i++)
//         {
//             if(intersection_map_[abscissa_idx][i].first == triangle_idx)
//                 skip = true;
//         }
//         if(! skip) {
//             storage.push_back(IntersectionLocal<1,2>(is));
//             intersection_map_[abscissa_idx].push_back(std::make_pair(
//                                                         triangle_idx,
//                                                         &(storage.back())
//                                                     ));
//             intersection_map_[triangle_idx].push_back(std::make_pair(
//                                                         abscissa_idx,
//                                                         &(storage.back())
//                                                     ));
//             DebugOut().fmt("1D-2D intersection [{} - {}]:\n",is.component_ele_idx(), is.bulk_ele_idx());
//             for(const IntersectionPointAux<1,2>& ip : is.points()) {
//                 //DebugOut() << ip;
//                 auto p = ip.coords(mesh->element(is.component_ele_idx()));
//                 DebugOut() << "[" << p[0] << " " << p[1] << " " << p[2] << "]\n";
//             }
//         }
//     }
//     END_TIMER("Intersection into storage");
}

void MixedMeshIntersections::compute_intersections_12_2(vector< IntersectionLocal< 1, 2 > >& storage)
{
    algorithm12_.compute_intersections_2(mesh->get_bih_tree());
//     DBGVAR(algorithm12_.intersectionaux_storage12_.size());
    
    START_TIMER("Intersection into storage");
    storage.reserve(algorithm12_.intersectionaux_storage12_.size());
    
    for(IntersectionAux<1,2> &is : algorithm12_.intersectionaux_storage12_) {
        store_intersection(storage, is);
//         DebugOut().fmt("1D-2D intersection [{} - {}]:\n",is.component_ele_idx(), is.bulk_ele_idx());
//         for(const IntersectionPointAux<1,2>& ip : is.points()) {
//             //DebugOut() << ip;
//             auto p = ip.coords(mesh->element(is.component_ele_idx()));
//             DebugOut() << "[" << p[0] << " " << p[1] << " " << p[2] << "]\n";
//         }
    }
    END_TIMER("Intersection into storage");
}

void MixedMeshIntersections::compute_intersections_12_1(vector< IntersectionLocal< 1, 2 > >& storage)
{
    algorithm12_.compute_intersections_1(mesh->get_bih_tree());
//     DBGVAR(algorithm12_.intersectionaux_storage12_.size());
    
    START_TIMER("Intersection into storage");
    storage.reserve(algorithm12_.intersectionaux_storage12_.size());
    
    for(IntersectionAux<1,2> &is : algorithm12_.intersectionaux_storage12_) {
        store_intersection(storage, is);
//         DebugOut().fmt("1D-2D intersection [{} - {}]:\n",is.component_ele_idx(), is.bulk_ele_idx());
//         for(const IntersectionPointAux<1,2>& ip : is.points()) {
//             //DebugOut() << ip;
//             auto p = ip.coords(mesh->element(is.component_ele_idx()));
//             DebugOut() << "[" << p[0] << " " << p[1] << " " << p[2] << "]\n";
//         }
    }
    END_TIMER("Intersection into storage");
}

void MixedMeshIntersections::compute_intersections(IntersectionType d)
{
    element_intersections_.resize(mesh->n_elements());
    
    // check whether the mesh is in plane only
    bool mesh_in_2d_only = false;
    auto bb = mesh->get_bih_tree().tree_box();
    for(uint axis = 0; axis < bb.dimension; axis++)
        if(bb.size(axis) < geometry_epsilon) mesh_in_2d_only = true;
    
    if(d & (IntersectionType::d13 | IntersectionType::d12_3)){
        START_TIMER("Intersections 1D-3D");
//         DebugOut() << "Intersection Algorithm d13\n";
        compute_intersections<1>(algorithm13_,intersection_storage13_);
        END_TIMER("Intersections 1D-3D");
    }
    append_to_index(intersection_storage13_);

    
    if(d & (IntersectionType::d23 | IntersectionType::d22 | IntersectionType::d12_3)){
        START_TIMER("Intersections 2D-3D");
//         DebugOut() << "Intersection Algorithm d23\n";
        compute_intersections<2>(algorithm23_,intersection_storage23_);
        END_TIMER("Intersections 2D-3D");
    }
    append_to_index(intersection_storage23_);
    

    if(d & IntersectionType::d22){
        START_TIMER("Intersections 2D-2D");
//         DebugOut() << "Intersection Algorithm d22\n";
        compute_intersections_22(intersection_storage22_);
        END_TIMER("Intersections 2D-2D");
    }

    if( mesh_in_2d_only){
        START_TIMER("Intersections 1D-2D (1)");
        if(d & IntersectionType::d12_1) compute_intersections_12_1(intersection_storage12_);
        END_TIMER("Intersections 1D-2D (1)");
    }
    // make sence only if some intersections in 3D are computed
    // TODO: this does NOT compute 1d-2d outside 3d bulk
    // NOTE: create input record in mesh to decide, whether compute also outside (means to call alg. 2)
    else if( ! intersection_storage13_.empty() &&
             ! intersection_storage23_.empty() &&
             (d & IntersectionType::d12_3)){
        START_TIMER("Intersections 1D-2D (3)");
        DebugOut() << "Intersection Algorithm d12_3\n";
        compute_intersections_12_3(intersection_storage12_);
        END_TIMER("Intersections 1D-2D (3)");
    }
    // otherwise compute 1d-2d in the most general case
    else if(d & IntersectionType::d12_2){
        START_TIMER("Intersections 1D-2D (2)");
        DebugOut() << "Intersection Algorithm d12_2\n";
        compute_intersections_12_2(intersection_storage12_);
        END_TIMER("Intersections 1D-2D (2)");
    }

    //ASSERT_EQ(intersection_storage13_.size(), 0);
    //ASSERT_EQ(intersection_storage23_.size(), 0);
    //ASSERT_EQ(intersection_storage22_.size(), 0);
    // compose master
    append_to_index(intersection_storage12_);

    // release temporary links from 3d elements
    for (auto elm : mesh->elements_range()) {
        if(elm->dim() == 3) element_intersections_[elm.idx()].clear();
    }



}



 
void MixedMeshIntersections::print_mesh_to_file_13(string name)
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
        for (auto nod : mesh->node_range()) {
            arma::vec3 _nod = *nod;
            fprintf(file,"%d %.16f %.16f %.16f\n", idx, _nod[0], _nod[1], _nod[2]);
            idx++;
        }

        for(unsigned int j = 0; j < intersection_storage13_.size();j++){
            IntersectionLocal<1,3> il = intersection_storage13_[j];
            ElementAccessor<3> el1D = mesh->element_accessor( il.component_ele_idx() );
//             ElementAccessor<3> el3D = mesh->element_accessor( il.bulk_ele_idx() );
            
            for(unsigned int k = 0; k < il.size();k++){
                number_of_nodes++;
                IntersectionPoint<1,3> IP13 = il[k];
                arma::vec3 global = IP13.coords(el1D);
                
//                 if(i == 0){
//                     _global = (IP13.local_bcoords_A())[0] * el1D.node(0)->point()
//                                                        +(IP13.local_bcoords_A())[1] * el1D.node(1)->point();
//                 }else{
//                     _global = (IP13.local_bcoords_B())[0] * el3D.node(0)->point()
//                                                        +(IP13.local_bcoords_B())[1] * el3D.node(1)->point()
//                                                        +(IP13.local_bcoords_B())[2] * el3D.node(2)->point()
//                                                        +(IP13.local_bcoords_B())[3] * el3D.node(3)->point();
//                 }

                fprintf(file,"%d %.16f %.16f %.16f\n", number_of_nodes, global[0], global[1], global[2]);
            }
        }

        fprintf(file,"$EndNodes\n");
        fprintf(file,"$Elements\n");
        fprintf(file,"%d\n", ((unsigned int)intersection_storage13_.size() + mesh->n_elements()) );

        for (auto elee : mesh->elements_range()) {
            if(elee->dim() == 3){
                int id1 = elee.node(0).idx() + 1;
                int id2 = elee.node(1).idx() + 1;
                int id3 = elee.node(2).idx() + 1;
                int id4 = elee.node(3).idx() + 1;

                fprintf(file,"%d 4 2 %d %d %d %d %d %d\n", elee.idx(), elee.region().id(), elee->pid(), id1, id2, id3, id4);
            }else if(elee->dim() == 2){
                int id1 = elee.node(0).idx() + 1;
                int id2 = elee.node(1).idx() + 1;
                int id3 = elee.node(2).idx() + 1;
                fprintf(file,"%d 2 2 %d %d %d %d %d\n", elee.idx(), elee.region().id(), elee->pid(), id1, id2, id3);

            }else if(elee->dim() == 1){
                int id1 = elee.node(0).idx() + 1;
                int id2 = elee.node(1).idx() + 1;
                fprintf(file,"%d 1 2 %d %d %d %d\n",elee.idx(), elee.region().id(), elee->pid(), id1, id2);
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

void MixedMeshIntersections::print_mesh_to_file_23(string name)
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
        for (auto nod : mesh->node_range()) {
            arma::vec3 _nod = *nod;
            fprintf(file,"%d %.16f %.16f %.16f\n", idx, _nod[0], _nod[1], _nod[2]);
            idx++;
        }

        for(unsigned int j = 0; j < intersection_storage23_.size();j++){
            
            IntersectionLocal<2,3> il = intersection_storage23_[j];
            ElementAccessor<3> el2D = mesh->element_accessor( il.component_ele_idx() );
//             ElementAccessor<3> el3D = mesh->element_accessor( il.bulk_ele_idx() );
            
            for(unsigned int k = 0; k < intersection_storage23_[j].size();k++){

                    number_of_nodes++;
                    IntersectionPoint<2,3> IP23 = il[k];
                    arma::vec3 global = IP23.coords(el2D);
//                     if(i == 0){
//                         _global = (IP23.local_bcoords_A())[0] * el2D.node(0)->point()
//                                                            +(IP23.local_bcoords_A())[1] * el2D.node(1)->point()
//                                                            +(IP23.local_bcoords_A())[2] * el2D.node(2)->point();
//                     }else{
//                         _global = (IP23.local_bcoords_B())[0] * el3D.node(0)->point()
//                                                            +(IP23.local_bcoords_B())[1] * el3D.node(1)->point()
//                                                            +(IP23.local_bcoords_B())[2] * el3D.node(2)->point()
//                                                            +(IP23.local_bcoords_B())[3] * el3D.node(3)->point();
//                     }
                    fprintf(file,"%d %.16f %.16f %.16f\n", number_of_nodes, global[0], global[1], global[2]);
            }
        }

        fprintf(file,"$EndNodes\n");
        fprintf(file,"$Elements\n");
        fprintf(file,"%d\n", (number_of_intersection_points + mesh->n_elements()) );

        for (auto elee : mesh->elements_range()) {
            if(elee->dim() == 3){
                int id1 = elee.node(0).idx() + 1;
                int id2 = elee.node(1).idx() + 1;
                int id3 = elee.node(2).idx() + 1;
                int id4 = elee.node(3).idx() + 1;

                fprintf(file,"%d 4 2 %d %d %d %d %d %d\n", elee.idx(), elee.region().id(), elee->pid(), id1, id2, id3, id4);
            }else if(elee->dim() == 2){
                int id1 = elee.node(0).idx() + 1;
                int id2 = elee.node(1).idx() + 1;
                int id3 = elee.node(2).idx() + 1;
                fprintf(file,"%d 2 2 %d %d %d %d %d\n", elee.idx(), elee.region().id(), elee->pid(), id1, id2, id3);

            }else{
                int id1 = elee.node(0).idx() + 1;
                int id2 = elee.node(1).idx() + 1;
                fprintf(file,"%d 1 2 %d %d %d %d\n",elee.idx(), elee.region().id(), elee->pid(), id1, id2);
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


