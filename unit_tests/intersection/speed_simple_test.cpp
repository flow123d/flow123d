/*
 * intersection_area_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: PE
 */
#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>
#include "system/global_defs.h"

#ifdef FLOW123D_RUN_UNIT_BENCHMARKS

#include <armadillo>

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "mesh/node_accessor.hh"
#include "mesh/range_wrapper.hh"
#include "io/msh_gmshreader.h"

#include "intersection/compute_intersection.hh"

#include "intersection/intersection_point_aux.hh"
#include "intersection/intersection_aux.hh"
#include "intersection/intersection_local.hh"

#include <dirent.h>

#include <ctime>
#include <cmath>

using namespace std;

static const unsigned int profiler_loop = 100;
static const unsigned int n_meshes = 10000;

// results - number of cases with number of ips 0-7
static unsigned int n_intersection[8] = {0, 0, 0, 0, 0, 0, 0, 0};
// results - number of pathologic cases with number of ips 0-7
// static unsigned int n_intersection_p[8] = {0, 0, 0, 0, 0, 0, 0, 0};

static void reset_statistics()
{
    for(unsigned int j=0; j<8; j++)
    {
        n_intersection[j] = 0;
//         n_intersection_p[j] = 0;
    } 
}

static void print_statistics()
{
    cout << "Results statistics:\nn_ips\tcount\tn_path\n-------------------------\n";
    for(unsigned int j=0; j<8; j++)
    {
        cout << setw(5) << j << "\t" << setw(5) << n_intersection[j] << endl;//"\t" << setw(6) << n_intersection_p[j] << endl;
    }
    cout << "-------------------------\n";
}

/// get random number in [0,1]
double unif(){
    return std::rand() / double(RAND_MAX);
}

/// get random number in [a,b]
double unif(double a, double b){
    return (b-a)*unif() + a;
}

/// set initial seed for rand using time
/// breaks the repeatibility of the test
void seed_rand(){
    std::srand(std::time(0));
}

// print testing mesh and save it to file of given name
void print_mesh(Mesh *mesh, string t_name = "random_mesh")
{
    unsigned int number_of_nodes = mesh->n_nodes();

    FILE * file;
    file = fopen((t_name.append(".msh")).c_str(),"w");
    
    fprintf(file, "$MeshFormat\n");
    fprintf(file, "2.2 0 8\n");
    fprintf(file, "$EndMeshFormat\n");
    fprintf(file, "$Nodes\n");
    fprintf(file, "%d\n", number_of_nodes);

    for (auto nod : mesh->node_range()) {
        arma::vec3 _nod = nod->point();
        fprintf(file,"%d %f %f %f\n", nod.idx()+1, _nod[0], _nod[1], _nod[2]);
    }

    fprintf(file,"$EndNodes\n");
    fprintf(file,"$Elements\n");
    fprintf(file,"%d\n",  mesh->n_elements() );

    for (auto elee : mesh->elements_range()) {
        if(elee->dim() == 3){
            int id1 = mesh->node_accessor(0).idx() + 1;
            int id2 = mesh->node_accessor(1).idx() + 1;
            int id3 = mesh->node_accessor(2).idx() + 1;
            int id4 = mesh->node_accessor(3).idx() + 1;

            fprintf(file,"%d 4 2 %d %d %d %d %d %d\n", elee.idx()+1, 3, elee->pid(), id1, id2, id3, id4);
        }else if(elee->dim() == 2){
            int id1 = mesh->node_accessor(0).idx() + 1;
            int id2 = mesh->node_accessor(1).idx() + 1;
            int id3 = mesh->node_accessor(2).idx() + 1;
            fprintf(file,"%d 2 2 %d %d %d %d %d\n", elee.idx()+1, 2, elee->pid(), id1, id2, id3);

        }else{
            int id1 = mesh->node_accessor(0).idx() + 1;
            int id2 = mesh->node_accessor(1).idx() + 1;
            fprintf(file,"%d 1 2 %d %d %d %d\n",elee.idx()+1, 1, elee->pid(), id1, id2);
        }
    }

    fprintf(file,"$EndElements\n");
    fclose(file);
}



// generates triangle vs tetrahedron mesh
template<unsigned int dimA, unsigned int dimB>
void generate_meshes(unsigned int N,
                     vector<Mesh*>& meshes)
{
    ASSERT(dimA <= dimB).error("Unsupported dimensions.");
    ASSERT(dimA != 3).error("Unsupported dimensions.");
    
    unsigned int nA = RefElement<dimA>::n_nodes, 
                 nB = RefElement<dimB>::n_nodes,
                 n_nodes = nA+nB;
//                  n_nodes = (nA+nB)*N;

    for (unsigned int i = 0; i < N; ++i) {
        Mesh* mesh = new Mesh();
        mesh->init_node_vector(n_nodes);
    	arma::vec3 point;
        for (unsigned int i = 0; i < n_nodes; ++i) {
            //generate random node
            point[0] = unif();
            point[1] = unif();
            point[2] = unif();
            mesh->add_node(i, point);
        }

        mesh->init_element_vector(2);

        std::vector<unsigned int> eleA_node_ids;
        for(unsigned int i =0; i < nA; i++)
        	eleA_node_ids.push_back(i);
        mesh->add_element(0, dimA, 1, 0, eleA_node_ids);

        std::vector<unsigned int> eleB_node_ids;
        for(unsigned int i =0; i < nB; i++)
        	eleB_node_ids.push_back(nA+i);
        // test tetrahedron node order
        if(dimB == 3)
        {
            double jac = arma::dot( arma::cross(mesh->node_accessor(nA+1)->point() - mesh->node_accessor(nA)->point(),
                                                mesh->node_accessor(nA+2)->point() - mesh->node_accessor(nA)->point()),
                                    mesh->node_accessor(nA+3)->point() - mesh->node_accessor(nA)->point());
            if( jac < 0)
            {
//                 DBGMSG("swap nodes: J = %f\n",jac);
                std::swap(eleB_node_ids[2], eleB_node_ids[3]);
            }
        }
        mesh->add_element(1, dimB, 2, 0, eleB_node_ids);
        
        mesh->side_nodes.resize(3);
        switch(dimA){
            case 1: mesh->side_nodes[0] = {{0},{1}}; break;
            case 2: mesh->side_nodes[1] = {{0,1},{0,2},{2,1}}; break;
            default: ASSERT(0)(dimA).error("Unsupported dimA");
        }
        switch(dimB){
            case 1: mesh->side_nodes[0] = {{nA+0},{nA+1}}; break;
            case 2: mesh->side_nodes[1] = {{nA+0, nA+1},{nA+0, nA+2},{nA+2, nA+1}}; break;
            case 3: mesh->side_nodes[2] = {{nA+0, nA+1, nA+2 }, { nA+0, nA+1, nA+3 }, { nA+0, nA+2, nA+3 }, { nA+1, nA+2, nA+3 }}; break;
            default: ASSERT(0)(dimA).error("Unsupported dimA");
        }
        
        meshes.push_back(mesh);
    }
}


template<unsigned int dimA, unsigned int dimB>
void compute_intersection(Mesh* mesh);

template<>
void compute_intersection<1,2>(Mesh* mesh)
{
	ElementAccessor<3> eleA = mesh->element_accessor(0);
	ElementAccessor<3> eleB = mesh->element_accessor(1);
    ASSERT_EQ(1, eleA->dim());
    ASSERT_EQ(2, eleB->dim());
    
    // compute intersection
    START_TIMER("Compute intersection");
   
    vector<Space<3>::Point> verticesA(2);
    vector<Space<3>::Point> verticesB(3);
    
    for(unsigned int i=0; i<2; i++) verticesA[i]=eleA.node(i)->point();
    for(unsigned int i=0; i<3; i++) verticesB[i]=eleB.node(i)->point();
     
    BoundingBox bbA(verticesA);
    BoundingBox bbB(verticesB);
    
    if(bbA.intersect(bbB)) {   
        START_TIMER("CI create");
        IntersectionAux<1,2> is(0, 1); //component_ele_idx, bulk_ele_idx
        ComputeIntersection<1,2> CI(eleA, eleB, mesh);
        END_TIMER("CI create");
        START_TIMER("CI compute");
        CI.compute_final(is.points());
        END_TIMER("CI compute");
        
        n_intersection[is.size()]++;
//         if(is.is_pathologic()) n_intersection_p[is.size()]++;
    }
    END_TIMER("Compute intersection");
}

template<unsigned int dimA, unsigned int dimB>
void compute_intersection(Mesh* mesh)
{
	ElementAccessor<3> eleA = mesh->element_accessor(0);
	ElementAccessor<3> eleB = mesh->element_accessor(1);
    ASSERT_EQ(dimA, eleA->dim());
    ASSERT_EQ(dimB, eleB->dim());
    // compute intersection
    START_TIMER("Compute intersection");
   
    vector<Space<3>::Point> verticesA(dimA+1);
    vector<Space<3>::Point> verticesB(dimB+1);
    
    for(unsigned int i=0; i<dimA+1; i++) verticesA[i]=eleA.node(i)->point();
    for(unsigned int i=0; i<dimB+1; i++) verticesB[i]=eleB.node(i)->point();
     
    BoundingBox bbA(verticesA);
    BoundingBox bbB(verticesB);
    
    if(bbA.intersect(bbB)) {   
        START_TIMER("CI create");
        IntersectionAux<dimA,dimB> is(0, 1); //component_ele_idx, bulk_ele_idx
        ComputeIntersection<dimA, dimB> CI(eleA, eleB, mesh);
        CI.init();
        END_TIMER("CI create");
        START_TIMER("CI compute");
        CI.compute(is);
        END_TIMER("CI compute");
        
        n_intersection[is.size()]++;
//         if(is.is_pathologic()) n_intersection_p[is.size()]++;
    }
    END_TIMER("Compute intersection");
}


// ***************************************************************************************************   1D-2D

TEST(speed_simple_12, all) {
    Profiler::instance();
    Profiler::instance()->set_task_info("Speed test for 1d-2d ComputeIntersection class.",2);
    
    reset_statistics();
    
    // create n random meshes triangle-tetrahedron in unit cube
    const unsigned int n = n_meshes;
    
    //seed_rand();
    vector<Mesh*> meshes;
    generate_meshes<1,2>(n,meshes);
    
    { START_TIMER("Speed test");
    // for each mesh, compute intersection area and compare with old NGH
    MessageOut() << "======== START ========\n";
    for(unsigned int i=0; i<n; i++)
    {       
            //MessageOut() << "================================================ %d\n",i);
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
            {
                compute_intersection<1,2>(meshes[i]);
            }
            //MessageOut() << "================================================\n";
    }
    MessageOut() << "======== FINISH ========\n";
    END_TIMER("Speed test"); }
    
    print_statistics();
    
    std::string profiler_file = "speed_simple_profiler_12.log";
    
    std::fstream fs;
    fs.open(profiler_file.c_str(), std::fstream::out);
    Profiler::instance()->output(PETSC_COMM_WORLD, fs);
    Profiler::uninitialize();
}

// ***************************************************************************************************   2D-2D

TEST(speed_simple_22, all) {
    Profiler::instance();
    Profiler::instance()->set_task_info("Speed test for 2d-2d ComputeIntersection class.",2);
    
    reset_statistics();
    
    // create n random meshes triangle-tetrahedron in unit cube
    const unsigned int n = n_meshes;
    
    //seed_rand();
    vector<Mesh*> meshes;
    generate_meshes<2,2>(n,meshes);
    
    { START_TIMER("Speed test");
    // for each mesh, compute intersection area and compare with old NGH
    MessageOut() << "======== START ========\n";
    for(unsigned int i=0; i<n; i++)
    {       
            //MessageOut() << "================================================ %d\n",i);
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
            {
                compute_intersection<2,2>(meshes[i]);
            }
            //MessageOut() << "================================================\n";
    }
    MessageOut() << "======== FINISH ========\n";
    END_TIMER("Speed test"); }
    
    print_statistics();
    
    std::string profiler_file = "speed_simple_profiler_22.log";
    
    std::fstream fs;
    fs.open(profiler_file.c_str(), std::fstream::out);
    Profiler::instance()->output(PETSC_COMM_WORLD, fs);
    Profiler::uninitialize();
}



// ***************************************************************************************************   1D-3D

TEST(speed_simple_13, all) {
    Profiler::instance();
    Profiler::instance()->set_task_info("Speed test for 1d-3d ComputeIntersection class.",2);
    
    reset_statistics();
    
    const unsigned int n = n_meshes;
    
    //seed_rand();
    vector<Mesh*> meshes;
    generate_meshes<1,3>(n,meshes);
    
    { START_TIMER("Speed test");
    // for each mesh, compute intersection area and compare with old NGH
    MessageOut() << "======== START ========\n";
    for(unsigned int i=0; i<n; i++)
    {       
            //MessageOut() << "================================================ %d\n",i);
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
            {
                compute_intersection<1,3>(meshes[i]);
            }
            //MessageOut() << "================================================\n";
    }
    MessageOut() << "======== FINISH ========\n";
    END_TIMER("Speed test"); }
    
    print_statistics();
    
    std::string profiler_file = "speed_simple_profiler_13.log";
    
    std::fstream fs;
    fs.open(profiler_file.c_str(), std::fstream::out);
    Profiler::instance()->output(PETSC_COMM_WORLD, fs);
    Profiler::uninitialize();
}



// ***************************************************************************************************   2D-3D

TEST(speed_simple_23, all) {
    Profiler::instance();
    Profiler::instance()->set_task_info("Speed test for 2d-3d ComputeIntersection class.",2);
    
    reset_statistics();
    
    // create n random meshes triangle-tetrahedron in unit cube
    const unsigned int n = n_meshes;
    
    //seed_rand();
    vector<Mesh*> meshes;
    generate_meshes<2,3>(n,meshes);
    
    { START_TIMER("Speed test");
    // for each mesh, compute intersection area and compare with old NGH
    MessageOut() << "======== START ========\n";
    for(unsigned int i=0; i<n; i++)
    {       
            //MessageOut() << "================================================ %d\n",i);
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
            {
                compute_intersection<2,3>(meshes[i]);
            }
            //MessageOut() << "================================================\n";
    }
    MessageOut() << "======== FINISH ========\n";
    END_TIMER("Speed test"); }
    
    print_statistics();
    
    std::string profiler_file = "speed_simple_profiler_23.log";
    
    std::fstream fs;
    fs.open(profiler_file.c_str(), std::fstream::out);
    Profiler::instance()->output(PETSC_COMM_WORLD, fs);
    Profiler::uninitialize();
}


#endif // FLOW123D_RUN_UNIT_BENCHMARKS
