#ifndef COMPUTE_INTERSECTION_TEST_H_
#define COMPUTE_INTERSECTION_TEST_H_

#include <dirent.h>

#include "system/file_path.hh"
#include "system/system.hh"
#include "mesh/mesh_types.hh"
#include "mesh/nodes.hh"
#include "mesh/elements.h"


static const std::vector<std::vector<unsigned int>> permutations_triangle = {
    {0,1,2},
    {1,0,2},
    {1,2,0},
    {0,2,1},
    {2,0,1},
    {2,1,0}};

static const std::vector<std::vector<unsigned int>> permutations_tetrahedron = {
    {0,1,2,3},
//         {0,1,3,2},  // the tab means permutation with negative jacobian
    {0,3,1,2},
//         {0,3,2,1},
    {0,2,3,1},
//         {0,2,1,3},
//         {1,0,2,3},
    {1,0,3,2},
//         {1,3,0,2},
    {1,3,2,0},
//         {1,2,3,0},
    {1,2,0,3},
//         {2,1,0,3},
    {2,1,3,0},
//         {2,3,1,0},
    {2,3,0,1},
//         {2,0,3,1},
    {2,0,1,3},
//         {3,1,2,0},
    {3,1,0,2},
//         {3,0,1,2},
    {3,0,2,1},
//         {3,2,0,1},
    {3,2,1,0}};

void read_files_from_dir(const string &dir_name,
                         const string &extension, 
                         std::vector<string> &filenames,
                         bool sort_filenames = true)
{
    // read mesh file names
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (dir_name.c_str())) != NULL) {
        // print all the files and directories within directory 
        MessageOut() << "Testing mesh files: \n";
        while ((ent = readdir (dir)) != NULL) {
            string fname = ent->d_name;
            // test extension ".msh"
            if(fname.size() >= 4)
            {
                string ext = fname.substr(fname.size()-3);
//                 xprintf(Msg,"%s\n",ext.c_str());
                if(ext == extension){
                    filenames.push_back(ent->d_name);
                    MessageOut() << ent->d_name << "\n";
                }
            }
        }
        closedir (dir);
    } else {
        ASSERT(0).error("Could not open directory with testing meshes.");
    }
    
    if(sort_filenames)
        std::sort(filenames.begin(), filenames.end(), less<string>());
}

void permute_tetrahedron(ElementFullIter ele, unsigned int p)
{
    ASSERT_DBG(ele->dim() == 3);
    ASSERT_DBG(p < permutations_tetrahedron.size());
    Node* tmp[4];
    for(unsigned int i=0; i<ele->n_nodes(); i++)
    {
        tmp[i] = ele->node[permutations_tetrahedron[p][i]];
    }
    for(unsigned int i=0; i<ele->n_nodes(); i++)
    {
        ele->node[i] = tmp[i];
//      ele->node[i]->point().print(cout);
    }
//  cout << p << ": jac = "  << ele->tetrahedron_jacobian() << endl;
}

void permute_triangle(ElementFullIter ele, unsigned int p)
{
    ASSERT_DBG(ele->dim() == 2);
    ASSERT_DBG(p < permutations_triangle.size());
    Node* tmp[3];
    for(unsigned int i=0; i<ele->n_nodes(); i++)
    {
        tmp[i] = ele->node[permutations_triangle[p][i]];
    }
    for(unsigned int i=0; i<ele->n_nodes(); i++)
    {
        ele->node[i] = tmp[i];
//      ele->node[i]->point().print(cout);
    }
//  cout << p << ": jac = "  << ele->tetrahedron_jacobian() << endl;
}


#endif // COMPUTE_INTERSECTION_TEST_H_
