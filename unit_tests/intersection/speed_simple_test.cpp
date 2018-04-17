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
#include "io/msh_gmshreader.h"

#include "mesh/ngh/include/point.h"
#include "mesh/ngh/include/intersection.h"

#include "intersection/simplex.hh"
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
static unsigned int n_intersection_p[8] = {0, 0, 0, 0, 0, 0, 0, 0};

static void reset_statistics()
{
    for(unsigned int j=0; j<8; j++)
    {
        n_intersection[j] = 0;
        n_intersection_p[j] = 0;
    } 
}

static void print_statistics()
{
    cout << "Results statistics:\nn_ips\tcount\tn_path\n-------------------------\n";
    for(unsigned int j=0; j<8; j++)
    {
        cout << setw(5) << j << "\t" << setw(5) << n_intersection[j] << "\t" << setw(6) << n_intersection_p[j] << endl;
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

    FOR_NODES(mesh, nod){
        arma::vec3 _nod = nod->point();
        fprintf(file,"%d %f %f %f\n", nod.id()+1, _nod[0], _nod[1], _nod[2]);
    }

    fprintf(file,"$EndNodes\n");
    fprintf(file,"$Elements\n");
    fprintf(file,"%d\n",  mesh->n_elements() );

    for (auto elee : mesh->bulk_elements_range()) {
        if(elee->dim() == 3){
            int id1 = mesh->node_vector.index(elee->node[0]) + 1;
            int id2 = mesh->node_vector.index(elee->node[1]) + 1;
            int id3 = mesh->node_vector.index(elee->node[2]) + 1;
            int id4 = mesh->node_vector.index(elee->node[3]) + 1;

            fprintf(file,"%d 4 2 %d %d %d %d %d %d\n", elee.idx()+1, 3, elee->pid(), id1, id2, id3, id4);
        }else if(elee->dim() == 2){
            int id1 = mesh->node_vector.index(elee->node[0]) + 1;
            int id2 = mesh->node_vector.index(elee->node[1]) + 1;
            int id3 = mesh->node_vector.index(elee->node[2]) + 1;
            fprintf(file,"%d 2 2 %d %d %d %d %d\n", elee.idx()+1, 2, elee->pid(), id1, id2, id3);

        }else{
            int id1 = mesh->node_vector.index(elee->node[0]) + 1;
            int id2 = mesh->node_vector.index(elee->node[1]) + 1;
            fprintf(file,"%d 1 2 %d %d %d %d\n",elee.idx()+1, 1, elee->pid(), id1, id2);
        }
    }

    fprintf(file,"$EndElements\n");
    fclose(file);
}



// generates triangle vs tetrahedron mesh
template<unsigned int dimA, unsigned int dimB>
void generate_meshes(unsigned int N,
                     vector<Simplex<dimA>>& eleA,
                     vector<Simplex<dimB>>& eleB,
                     vector<Space<3>::Point*>& nodes)
{
    ASSERT(dimA <= dimB).error("Unsupported dimensions.");
    ASSERT(dimA != 3).error("Unsupported dimensions.");
    
    unsigned int nA = RefElement<dimA>::n_nodes, 
                 nB = RefElement<dimB>::n_nodes,
                 n_nodes = (nA+nB)*N;

    eleA.resize(N);
    eleB.resize(N);
    nodes.reserve(n_nodes);
    
    //generate random nodes
    for (unsigned int i = 0; i < n_nodes; ++i) {
        //generate random node
        nodes.push_back( new Space<3>::Point({unif(), unif(), unif()}) );
    }

    // set nodes to simplices: [eleA nodes; eleB nodes]
    for (unsigned int i = 0; i < N; ++i) {
        unsigned int offsetA = i*nA;           //first part of nodes
        unsigned int offsetB = N*nA + i*nB;    //second part of nodes
        
        eleA[i].set_simplices(nodes.data() + offsetA);
        eleB[i].set_simplices(nodes.data() + offsetB);
        
        if(dimB == 3) {
            double jac = arma::dot( arma::cross(*nodes[offsetB+1] - *nodes[offsetB], 
                                                *nodes[offsetB+2] - *nodes[offsetB]),
                                    *nodes[offsetB+3] - *nodes[offsetB]);
            if( jac < 0)
            {
//                 DBGMSG("swap nodes: J = %f\n",jac);
                std::swap(*nodes[offsetB+2], *nodes[offsetB+3]);
            }
//             DBGMSG("J = %f\n",jac);
        }
    }
}


template<unsigned int dimA, unsigned int dimB>
void compute_intersection(Simplex<dimA>& eleA,
                          Simplex<dimB>& eleB);

template<>
void compute_intersection<1,2>(Simplex<1>& eleA,
                               Simplex<2>& eleB)
{
    // compute intersection
    
    START_TIMER("Compute intersection");
   
    vector<Space<3>::Point> verticesA(2);
    vector<Space<3>::Point> verticesB(3);
    
    for(unsigned int i=0; i<2; i++) verticesA[i]=eleA.node(i).point_coordinates();
    for(unsigned int i=0; i<3; i++) verticesB[i]=eleB.node(i).point_coordinates();
     
    BoundingBox bbA(verticesA);
    BoundingBox bbB(verticesB);
    
    if(bbA.intersect(bbB)) {   
        START_TIMER("CI create");
        IntersectionAux<1,2> is(0, 1, 0); //component_ele_idx, bulk_ele_idx, component_idx
        ComputeIntersection<Simplex<1>, Simplex<2>> CI(eleA, eleB, mesh);
        END_TIMER("CI create");
        START_TIMER("CI compute");
        CI.compute_final(is.points());
        END_TIMER("CI compute");
        
        n_intersection[is.size()]++;
        if(is.is_pathologic()) n_intersection_p[is.size()]++;
    }
    END_TIMER("Compute intersection");
}

template<unsigned int dimA, unsigned int dimB>
void compute_intersection(Simplex<dimA>& eleA,
                          Simplex<dimB>& eleB)
{
    // compute intersection
    START_TIMER("Compute intersection");
   
    vector<Space<3>::Point> verticesA(dimA+1);
    vector<Space<3>::Point> verticesB(dimB+1);
    
    for(unsigned int i=0; i<dimA+1; i++) verticesA[i]=eleA.node(i).point_coordinates();
    for(unsigned int i=0; i<dimB+1; i++) verticesB[i]=eleB.node(i).point_coordinates();
     
    BoundingBox bbA(verticesA);
    BoundingBox bbB(verticesB);
    
    if(bbA.intersect(bbB)) {   
        START_TIMER("CI create");
        IntersectionAux<dimA,dimB> is(0, 1, 0); //component_ele_idx, bulk_ele_idx, component_idx
        ComputeIntersection<Simplex<dimA>, Simplex<dimB>> CI(eleA, eleB, mesh);
        CI.init();
        END_TIMER("CI create");
        START_TIMER("CI compute");
        CI.compute(is);
        END_TIMER("CI compute");
        
        n_intersection[is.size()]++;
        if(is.is_pathologic()) n_intersection_p[is.size()]++;
    }
    END_TIMER("Compute intersection");
}


void compute_intersection_ngh_12(Simplex<1>& eleA,
                                 Simplex<2>& eleB)
{
    TAbscissa tabs;
    TTriangle ttr;
    IntersectionLocal* is;

    START_TIMER("Compute intersection NGH");

    tabs.SetPoints(TPoint(eleA.node(0).point_coordinates()(0), eleA.node(0).point_coordinates()(1), eleA.node(0).point_coordinates()(2)),
                   TPoint(eleA.node(1).point_coordinates()(0), eleA.node(1).point_coordinates()(1), eleA.node(1).point_coordinates()(2)));

    ttr.SetPoints(TPoint(eleB.node(0).point_coordinates()(0), eleB.node(0).point_coordinates()(1), eleB.node(0).point_coordinates()(2)),
                  TPoint(eleB.node(1).point_coordinates()(0), eleB.node(1).point_coordinates()(1), eleB.node(1).point_coordinates()(2)),
                  TPoint(eleB.node(2).point_coordinates()(0), eleB.node(2).point_coordinates()(1), eleB.node(2).point_coordinates()(2)) );

    GetIntersection(tabs, ttr, is);
    END_TIMER("Compute intersection NGH");
}

void compute_intersection_ngh_22(Simplex<2>& eleA,
                                 Simplex<2>& eleB)
{
    TTriangle ttrA, ttrB;
    TIntersectionType it = unknown;
    double area;

    START_TIMER("Compute intersection NGH");

    ttrA.SetPoints(TPoint(eleA.node(0).point_coordinates()(0), eleA.node(0).point_coordinates()(1), eleA.node(0).point_coordinates()(2)),
                   TPoint(eleA.node(1).point_coordinates()(0), eleA.node(1).point_coordinates()(1), eleA.node(1).point_coordinates()(2)),
                   TPoint(eleA.node(2).point_coordinates()(0), eleA.node(2).point_coordinates()(1), eleA.node(2).point_coordinates()(2)) );

    ttrB.SetPoints(TPoint(eleB.node(0).point_coordinates()(0), eleB.node(0).point_coordinates()(1), eleB.node(0).point_coordinates()(2)),
                   TPoint(eleB.node(1).point_coordinates()(0), eleB.node(1).point_coordinates()(1), eleB.node(1).point_coordinates()(2)),
                   TPoint(eleB.node(2).point_coordinates()(0), eleB.node(2).point_coordinates()(1), eleB.node(2).point_coordinates()(2)) );

    GetIntersection(ttrA, ttrB, it, area);
    END_TIMER("Compute intersection NGH");
} 

void compute_intersection_ngh_13(Simplex<1>& eleA,
                                 Simplex<3>& eleB)
{
    double length;
    
    { START_TIMER("Compute intersection NGH");
        
    TPoint p1 = TPoint(eleA.node(0).point_coordinates()(0), eleA.node(0).point_coordinates()(1), eleA.node(0).point_coordinates()(2)),
           p2 = TPoint(eleA.node(1).point_coordinates()(0), eleA.node(1).point_coordinates()(1), eleA.node(1).point_coordinates()(2)),
           p3 = TPoint(eleB.node(0).point_coordinates()(0), eleB.node(0).point_coordinates()(1), eleB.node(0).point_coordinates()(2)),
           p4 = TPoint(eleB.node(1).point_coordinates()(0), eleB.node(1).point_coordinates()(1), eleB.node(1).point_coordinates()(2)),
           p5 = TPoint(eleB.node(2).point_coordinates()(0), eleB.node(2).point_coordinates()(1), eleB.node(2).point_coordinates()(2)),
           p6 = TPoint(eleB.node(3).point_coordinates()(0), eleB.node(3).point_coordinates()(1), eleB.node(3).point_coordinates()(2));

    TAbscissa tabs(p1,p2);
    TTetrahedron tte(p3, p4, p5, p6);
    TIntersectionType it = Intersections::line;
    
    GetIntersection(tabs, tte, it, length);
    END_TIMER("Compute intersection NGH"); }
} 

void compute_intersection_ngh_23(Simplex<2>& eleA,
                                 Simplex<3>& eleB)
{
    double area;
    
    {START_TIMER("Compute intersection NGH");
        
    TPoint p1 = TPoint(eleA.node(0).point_coordinates()(0), eleA.node(0).point_coordinates()(1), eleA.node(0).point_coordinates()(2)),
           p2 = TPoint(eleA.node(1).point_coordinates()(0), eleA.node(1).point_coordinates()(1), eleA.node(1).point_coordinates()(2)),
           p3 = TPoint(eleA.node(2).point_coordinates()(0), eleA.node(2).point_coordinates()(1), eleA.node(2).point_coordinates()(2)),
           p4 = TPoint(eleB.node(0).point_coordinates()(0), eleB.node(0).point_coordinates()(1), eleB.node(0).point_coordinates()(2)),
           p5 = TPoint(eleB.node(1).point_coordinates()(0), eleB.node(1).point_coordinates()(1), eleB.node(1).point_coordinates()(2)),
           p6 = TPoint(eleB.node(2).point_coordinates()(0), eleB.node(2).point_coordinates()(1), eleB.node(2).point_coordinates()(2)),
           p7 = TPoint(eleB.node(3).point_coordinates()(0), eleB.node(3).point_coordinates()(1), eleB.node(3).point_coordinates()(2));
           
        
    TTriangle ttr(p1,p2,p3);
    TTetrahedron tte(p4,p5,p6,p7);
    TIntersectionType it = Intersections::area;

    GetIntersection(ttr, tte, it, area);
    END_TIMER("Compute intersection NGH");}
} 




// ***************************************************************************************************   1D-2D

TEST(speed_simple_12, all) {
    Profiler::initialize();
    Profiler::instance()->set_task_info("Speed test for 1d-2d ComputeIntersection class.",2);
    
    reset_statistics();
    
    // create n random meshes triangle-tetrahedron in unit cube
    const unsigned int n = n_meshes;
    
    //seed_rand();
    vector<Simplex<1>> eleA;
    vector<Simplex<2>> eleB;
    vector<Space<3>::Point*> nodes;
    generate_meshes<1,2>(n,eleA, eleB, nodes);
    
    
    { START_TIMER("Speed test NGH");
    // for each mesh, compute intersection area and compare with old NGH
    MessageOut() << "======== NGH ========\n";
    for(unsigned int i=0; i<n; i++)
    {       
            //MessageOut() << "================================================ %d\n",i);
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
            {
                compute_intersection_ngh_12(eleA[i], eleB[i]);
            }
            //MessageOut() << "================================================\n";
    }
    MessageOut() << "======== NGH end ========\n";
    END_TIMER("Speed test NGH"); }
    
    { START_TIMER("Speed test");
    // for each mesh, compute intersection area and compare with old NGH
    MessageOut() << "======== NEW ========\n";
    for(unsigned int i=0; i<n; i++)
    {       
            //MessageOut() << "================================================ %d\n",i);
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
            {
                compute_intersection<1,2>(eleA[i], eleB[i]);
            }
            //MessageOut() << "================================================\n";
    }
    MessageOut() << "======== NEW end ========\n";
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
    Profiler::initialize();
    Profiler::instance()->set_task_info("Speed test for 2d-2d ComputeIntersection class.",2);
    
    reset_statistics();
    
    // create n random meshes triangle-tetrahedron in unit cube
    const unsigned int n = n_meshes;
    
    //seed_rand();
    vector<Simplex<2>> eleA;
    vector<Simplex<2>> eleB;
    vector<Space<3>::Point*> nodes;
    generate_meshes<2,2>(n,eleA, eleB, nodes);
    
    
    { START_TIMER("Speed test NGH");
    // for each mesh, compute intersection area and compare with old NGH
    MessageOut() << "======== NGH ========\n";
    for(unsigned int i=0; i<n; i++)
    {       
            //MessageOut() << "================================================ %d\n",i);
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
            {
                compute_intersection_ngh_22(eleA[i], eleB[i]);
            }
            //MessageOut() << "================================================\n";
    }
    MessageOut() << "======== NGH end ========\n";
    END_TIMER("Speed test NGH"); }
    
    { START_TIMER("Speed test");
    // for each mesh, compute intersection area and compare with old NGH
    MessageOut() << "======== NEW ========\n";
    for(unsigned int i=0; i<n; i++)
    {       
            //MessageOut() << "================================================ %d\n",i);
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
            {
                compute_intersection<2,2>(eleA[i], eleB[i]);
            }
            //MessageOut() << "================================================\n";
    }
    MessageOut() << "======== NEW end ========\n";
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
    Profiler::initialize();
    Profiler::instance()->set_task_info("Speed test for 1d-3d ComputeIntersection class.",2);
    
    reset_statistics();
    
    const unsigned int n = n_meshes;
    
    //seed_rand();
    vector<Simplex<1>> eleA;
    vector<Simplex<3>> eleB;
    vector<Space<3>::Point*> nodes;
    generate_meshes<1,3>(n,eleA, eleB, nodes);

    
    { START_TIMER("Speed test NGH");
    // for each mesh, compute intersection area and compare with old NGH
    MessageOut() << "======== NGH ========\n";
    for(unsigned int i=0; i<n; i++)
    {       
            //MessageOut() << "================================================ %d\n",i);
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
            {
                compute_intersection_ngh_13(eleA[i], eleB[i]);
            }
            //MessageOut() << "================================================\n";
    }
    MessageOut() << "======== NGH end ========\n";
    END_TIMER("Speed test NGH"); }
    
    { START_TIMER("Speed test");
    // for each mesh, compute intersection area and compare with old NGH
    MessageOut() << "======== NEW ========\n";
    for(unsigned int i=0; i<n; i++)
    {       
            //MessageOut() << "================================================ %d\n",i);
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
            {
                compute_intersection(eleA[i], eleB[i]);
            }
            //MessageOut() << "================================================\n";
    }
    MessageOut() << "======== NEW end ========\n";
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
    Profiler::initialize();
    Profiler::instance()->set_task_info("Speed test for 2d-3d ComputeIntersection class.",2);
    
    reset_statistics();
    
    // create n random meshes triangle-tetrahedron in unit cube
    const unsigned int n = n_meshes;
    
    //seed_rand();
    vector<Simplex<2>> eleA;
    vector<Simplex<3>> eleB;
    vector<Space<3>::Point*> nodes;
    generate_meshes<2,3>(n,eleA, eleB, nodes);
    
    
    { START_TIMER("Speed test NGH");
    // for each mesh, compute intersection area and compare with old NGH
    MessageOut() << "======== NGH ========\n";
    for(unsigned int i=0; i<n; i++)
    {       
            //MessageOut() << "================================================ %d\n",i);
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
            {
                compute_intersection_ngh_23(eleA[i], eleB[i]);
            }
            //MessageOut() << "================================================\n";
    }
    MessageOut() << "======== NGH end ========\n";
    END_TIMER("Speed test NGH"); }
    
    { START_TIMER("Speed test");
    // for each mesh, compute intersection area and compare with old NGH
    MessageOut() << "======== NEW ========\n";
    for(unsigned int i=0; i<n; i++)
    {       
            //MessageOut() << "================================================ %d\n",i);
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
            {
                compute_intersection<2,3>(eleA[i], eleB[i]);
            }
            //MessageOut() << "================================================\n";
    }
    MessageOut() << "======== NEW end ========\n";
    END_TIMER("Speed test"); }
    
    print_statistics();
    
    std::string profiler_file = "speed_simple_profiler_23.log";
    
    std::fstream fs;
    fs.open(profiler_file.c_str(), std::fstream::out);
    Profiler::instance()->output(PETSC_COMM_WORLD, fs);
    Profiler::uninitialize();
}


#endif // FLOW123D_RUN_UNIT_BENCHMARKS
