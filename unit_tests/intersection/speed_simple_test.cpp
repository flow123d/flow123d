/*
 * intersection_area_test.cpp
 *
 *  Created on: Nov 14, 2014
 *      Author: PE
 */
#define TEST_USE_PETSC
#include <flow_gtest_mpi.hh>

#include <armadillo>

#include "system/system.hh"
#include "system/sys_profiler.hh"
#include "system/file_path.hh"
#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"

#include "mesh/ngh/include/point.h"
#include "mesh/ngh/include/intersection.h"

#include "intersection/simplex.h"
#include "intersection/computeintersection.h"

#include "intersection/intersectionpoint.h"
#include "intersection/intersectionaux.h"
#include "intersection/intersection_local.h"

#include <dirent.h>

#include <ctime>
#include <cmath>

using namespace std;

static const unsigned int profiler_loop = 100;
static const unsigned int n_meshes = 100000;

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

    FOR_ELEMENTS(mesh, elee){
        if(elee->dim() == 3){
            int id1 = mesh->node_vector.index(elee->node[0]) + 1;
            int id2 = mesh->node_vector.index(elee->node[1]) + 1;
            int id3 = mesh->node_vector.index(elee->node[2]) + 1;
            int id4 = mesh->node_vector.index(elee->node[3]) + 1;

            fprintf(file,"%d 4 2 %d %d %d %d %d %d\n", elee.id()+1, 3, elee->pid, id1, id2, id3, id4);
        }else if(elee->dim() == 2){
            int id1 = mesh->node_vector.index(elee->node[0]) + 1;
            int id2 = mesh->node_vector.index(elee->node[1]) + 1;
            int id3 = mesh->node_vector.index(elee->node[2]) + 1;
            fprintf(file,"%d 2 2 %d %d %d %d %d\n", elee.id()+1, 2, elee->pid, id1, id2, id3);

        }else{
            int id1 = mesh->node_vector.index(elee->node[0]) + 1;
            int id2 = mesh->node_vector.index(elee->node[1]) + 1;
            fprintf(file,"%d 1 2 %d %d %d %d\n",elee.id()+1, 1, elee->pid, id1, id2);
        }
    }

    fprintf(file,"$EndElements\n");
    fclose(file);
}



// generates triangle vs tetrahedron mesh
template<unsigned int dimA, unsigned int dimB>
void generate_mesh(Mesh *mesh)
{
    ASSERT(dimA <= dimB, "Unsupported dimensions.");
    ASSERT(dimA != 3, "Unsupported dimensions.");
    
    mesh->node_vector.reserve(10);
    unsigned int nA = RefElement<dimA>::n_nodes, 
                 nB = RefElement<dimB>::n_nodes;
    
    for (unsigned int i = 0; i < nA+nB; ++i) {
        NodeFullIter node = mesh->node_vector.add_item(i);

        node->point()(0) = unif();
        node->point()(1) = unif();
        node->point()(2) = unif();
    }
    
    mesh->element.reserve(2);
    unsigned int ele_id = 0;

    // element A
    Element *ele = mesh->element.add_item(ele_id);
    RegionIdx reg_idx;
    ele->init(dimA, mesh, reg_idx);
    ele->pid = 0;

    for (unsigned int ni = 0; ni < nA; ++ni) {
        unsigned int node_id = ni;
        NodeIter node = mesh->node_vector.find_id( node_id );
        INPUT_CHECK( node!=mesh->node_vector.end() ,
                "Unknown node id %d in specification of element.\n",
                node_id);
        ele->node[ni] = node;
    }
    
    // element B
    ele_id++;
    ele = mesh->element.add_item(ele_id);
    ele->init(dimB, mesh, reg_idx);
    ele->pid = 0;

    for (unsigned int ni = 0; ni < nB; ++ni) {
        unsigned int node_id = nA + ni;
        NodeIter node = mesh->node_vector.find_id( node_id );
        INPUT_CHECK( node!=mesh->node_vector.end() ,
                "Unknown node id %d in specification of element.\n",
                node_id);
        ele->node[ni] = node;
    }
    
    if(dimB == 3) {
        double jac = ele->tetrahedron_jacobian();
        if( jac < 0)
        {
            DBGMSG("swap nodes: J = %f\n",jac);
            std::swap(*ele->node[2], *ele->node[3]);
        }
        DBGMSG("J = %f\n",ele->tetrahedron_jacobian());
    }
    
    mesh->setup_topology();
}


/// Auxiliary function that translates @p ElementFullIter to @p Simplex<simplex_dim>.
template<unsigned int simplex_dim>
void update_simplex(const ElementFullIter &element, computeintersection::Simplex<simplex_dim> & simplex)
{
    arma::vec3 *field_of_points[simplex_dim+1];
    for(unsigned int i=0; i < simplex_dim+1; i++)
        field_of_points[i]= &(element->node[i]->point());
    simplex.set_simplices(field_of_points);
}

template<unsigned int dimA, unsigned int dimB>
void compute_intersection(Mesh *mesh)
{
    // compute intersection
    std::vector< unsigned int > prolongation_table;
    
    ElementFullIter elmA = mesh->element(0);
    ElementFullIter elmB = mesh->element(1);
    
    START_TIMER("Compute intersection");
   
    BoundingBox bbA = elmA->bounding_box();
    BoundingBox bbB = elmB->bounding_box();
    
    if(bbA.intersect(bbB)) {
    
        computeintersection::Simplex<dimA> comp_simplex;
        computeintersection::Simplex<dimB> bulk_simplex;
        update_simplex(elmA, comp_simplex);
        update_simplex(elmB, bulk_simplex);
    
        START_TIMER("CI create");
        computeintersection::IntersectionAux<dimA,dimB> is(0, 1, 0); //component_ele_idx, bulk_ele_idx, component_idx
        computeintersection::ComputeIntersection<computeintersection::Simplex<dimA>, computeintersection::Simplex<dimB>> CI(comp_simplex, bulk_simplex);
        CI.init();
        END_TIMER("CI create");
        START_TIMER("CI compute");
        CI.compute(is, prolongation_table);
        END_TIMER("CI compute");
        
        n_intersection[is.size()]++;
        if(is.is_pathologic()) n_intersection_p[is.size()]++;
    }
    END_TIMER("Compute intersection");
}


void compute_intersection_ngh_12(Mesh *mesh)
{
    TAbscissa tabs;
    TTriangle ttr;
    IntersectionLocal* is;
    
    ElementFullIter elmA = mesh->element(0);
    ElementFullIter elmB = mesh->element(1);

    START_TIMER("Compute intersection NGH");

    tabs.SetPoints(TPoint(elmA->node[0]->point()(0), elmA->node[0]->point()(1), elmA->node[0]->point()(2)),
                   TPoint(elmA->node[1]->point()(0), elmA->node[1]->point()(1), elmA->node[1]->point()(2)));

    ttr.SetPoints(TPoint(elmB->node[0]->point()(0), elmB->node[0]->point()(1), elmB->node[0]->point()(2)),
                  TPoint(elmB->node[1]->point()(0), elmB->node[1]->point()(1), elmB->node[1]->point()(2)),
                  TPoint(elmB->node[2]->point()(0), elmB->node[2]->point()(1), elmB->node[2]->point()(2)) );

    GetIntersection(tabs, ttr, is);
    END_TIMER("Compute intersection NGH");
}

void compute_intersection_ngh_22(Mesh *mesh)
{
    TTriangle ttrA, ttrB;
    TIntersectionType it = unknown;
    double area;
    
    ElementFullIter elmA = mesh->element(0);
    ElementFullIter elmB = mesh->element(1);

    START_TIMER("Compute intersection NGH");

    ttrA.SetPoints(TPoint(elmA->node[0]->point()(0), elmA->node[0]->point()(1), elmA->node[0]->point()(2)),
                   TPoint(elmA->node[1]->point()(0), elmA->node[1]->point()(1), elmA->node[1]->point()(2)),
                   TPoint(elmA->node[2]->point()(0), elmA->node[2]->point()(1), elmA->node[2]->point()(2)) );

    ttrB.SetPoints(TPoint(elmB->node[0]->point()(0), elmB->node[0]->point()(1), elmB->node[0]->point()(2)),
                   TPoint(elmB->node[1]->point()(0), elmB->node[1]->point()(1), elmB->node[1]->point()(2)),
                   TPoint(elmB->node[2]->point()(0), elmB->node[2]->point()(1), elmB->node[2]->point()(2)) );

    GetIntersection(ttrA, ttrB, it, area);
    END_TIMER("Compute intersection NGH");
} 

void compute_intersection_ngh_13(Mesh *mesh)
{
    double length;
    ElementFullIter elmA = mesh->element(0);
    ElementFullIter elmB = mesh->element(1);
    
    { START_TIMER("Compute intersection NGH");
        
    TPoint p1 = TPoint(elmA->node[0]->point()(0), elmA->node[0]->point()(1), elmA->node[0]->point()(2)),
           p2 = TPoint(elmA->node[1]->point()(0), elmA->node[1]->point()(1), elmA->node[1]->point()(2)),
           p3 = TPoint(elmB->node[0]->point()(0), elmB->node[0]->point()(1), elmB->node[0]->point()(2)),
           p4 = TPoint(elmB->node[1]->point()(0), elmB->node[1]->point()(1), elmB->node[1]->point()(2)),
           p5 = TPoint(elmB->node[2]->point()(0), elmB->node[2]->point()(1), elmB->node[2]->point()(2)),
           p6 = TPoint(elmB->node[3]->point()(0), elmB->node[3]->point()(1), elmB->node[3]->point()(2));

    TAbscissa tabs(p1,p2);
    TTetrahedron tte(p3, p4, p5, p6);
    TIntersectionType it = Intersections::line;
    
    GetIntersection(tabs, tte, it, length);
    END_TIMER("Compute intersection NGH"); }
} 

void compute_intersection_ngh_23(Mesh *mesh)
{
    double area;
    ElementFullIter elmA = mesh->element(0);
    ElementFullIter elmB = mesh->element(1);
    
    {START_TIMER("Compute intersection NGH");
        
    TPoint p1 = TPoint(elmA->node[0]->point()(0), elmA->node[0]->point()(1), elmA->node[0]->point()(2)),
           p2 = TPoint(elmA->node[1]->point()(0), elmA->node[1]->point()(1), elmA->node[1]->point()(2)),
           p3 = TPoint(elmA->node[2]->point()(0), elmA->node[2]->point()(1), elmA->node[2]->point()(2)),
           p4 = TPoint(elmB->node[0]->point()(0), elmB->node[0]->point()(1), elmB->node[0]->point()(2)),
           p5 = TPoint(elmB->node[1]->point()(0), elmB->node[1]->point()(1), elmB->node[1]->point()(2)),
           p6 = TPoint(elmB->node[2]->point()(0), elmB->node[2]->point()(1), elmB->node[2]->point()(2)),
           p7 = TPoint(elmB->node[3]->point()(0), elmB->node[3]->point()(1), elmB->node[3]->point()(2));
           
        
    TTriangle ttr(p1,p2,p3);
    TTetrahedron tte(p4,p5,p6,p7);
    TIntersectionType it = Intersections::area;

    GetIntersection(ttr, tte, it, area);
    END_TIMER("Compute intersection NGH");}
} 






// ***************************************************************************************************   1D-3D

TEST(speed_simple_13, all) {
    Profiler::initialize();
    Profiler::instance()->set_task_info("Speed test for 1d-3d ComputeIntersection class.",2);
    
    reset_statistics();
    
    // create n random meshes triangle-tetrahedron in unit cube
    const unsigned int n = n_meshes;
    Mesh meshes[n];
    
    //seed_rand();
    for(unsigned int i=0; i<n; i++)
    {
        generate_mesh<1,3>(&meshes[i]);
        
//         stringstream stream;
//         stream << "random_meshes/rand_mesh_" << i << ".msh";
//         print_mesh(&meshes[i], stream.str());
    }
    
    { START_TIMER("Speed test NGH");
    // for each mesh, compute intersection area and compare with old NGH
    xprintf(Msg, "======== NGH ========\n");
    for(unsigned int i=0; i<n; i++)
    {       
            //xprintf(Msg, "================================================ %d\n",i);
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
            {
                compute_intersection_ngh_13(&meshes[i]);
            }
            //xprintf(Msg, "================================================\n");
    }
    xprintf(Msg, "======== NGH end ========\n");
    END_TIMER("Speed test NGH"); }
    
    { START_TIMER("Speed test");
    // for each mesh, compute intersection area and compare with old NGH
    xprintf(Msg, "======== NEW ========\n");
    for(unsigned int i=0; i<n; i++)
    {       
            //xprintf(Msg, "================================================ %d\n",i);
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
            {
                compute_intersection<1,3>(&meshes[i]);
            }
            //xprintf(Msg, "================================================\n");
    }
    xprintf(Msg, "======== NEW end ========\n");
    END_TIMER("Speed test"); }
    
    print_statistics();
    
    std::string profiler_file = "speed_simple_profiler_13.log";
    
    std::fstream fs;
    fs.open(profiler_file.c_str(), std::fstream::out | std::fstream::app);
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
    Mesh meshes[n];
    
    //seed_rand();
    for(unsigned int i=0; i<n; i++)
    {
        generate_mesh<2,3>(&meshes[i]);
        
//         stringstream stream;
//         stream << "random_meshes/rand_mesh_" << i << ".msh";
//         print_mesh(&meshes[i], stream.str());
    }
    
    { START_TIMER("Speed test NGH");
    // for each mesh, compute intersection area and compare with old NGH
    xprintf(Msg, "======== NGH ========\n");
    for(unsigned int i=0; i<n; i++)
    {       
            //xprintf(Msg, "================================================ %d\n",i);
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
            {
                compute_intersection_ngh_23(&meshes[i]);
            }
            //xprintf(Msg, "================================================\n");
    }
    xprintf(Msg, "======== NGH end ========\n");
    END_TIMER("Speed test NGH"); }
    
    { START_TIMER("Speed test");
    // for each mesh, compute intersection area and compare with old NGH
    xprintf(Msg, "======== NEW ========\n");
    for(unsigned int i=0; i<n; i++)
    {       
            //xprintf(Msg, "================================================ %d\n",i);
            for(unsigned int loop = 0; loop < profiler_loop; loop++)
            {
                compute_intersection<2,3>(&meshes[i]);
            }
            //xprintf(Msg, "================================================\n");
    }
    xprintf(Msg, "======== NEW end ========\n");
    END_TIMER("Speed test"); }
    
    print_statistics();
    
    std::string profiler_file = "speed_simple_profiler_23.log";
    
    std::fstream fs;
    fs.open(profiler_file.c_str(), std::fstream::out | std::fstream::app);
    Profiler::instance()->output(PETSC_COMM_WORLD, fs);
    Profiler::uninitialize();
}