/*!
 *
ï»¿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 *
 * @file    dfs_topo_sort.hh
 * @brief
 */

#ifndef DFS_TOPO_SORT_HH_
#define DFS_TOPO_SORT_HH_

#include <list>
#include <stack>
#include <vector>

using namespace std;

/**
 * Class to represent a graph and allows to perform sorting.
 *
 * DFS topological sorting is used.
 */
class DfsTopoSort {
    /// Number of vertices
	unsigned int nv_;

    /// Pointer to an array containing adjacency lists
    list<unsigned int>* adj_;

    /// A function used by topological_sort
    void topological_sort_util(unsigned int v, bool visited[],
                             stack<unsigned int>& stc);

public:
    /// Constructor
    DfsTopoSort(unsigned int nv)
    : nv_(nv)
    {
    	adj_ = new list<unsigned int>[nv];
    }

    /// Destructor
    ~DfsTopoSort()
    {
        delete[] adj_;
    }

    /// Clear content of adjacency lists
    inline void clear()
    {
        for (uint i=0; i<nv_; ++i)
    	    adj_[i].clear();
    }

    /**
     * Add an edge to graph
     *
     * Vertex 'w' is descendant of 'v'
     */
    inline void add_edge(unsigned int v, unsigned int w)
    {
        // Add w to v’s list.
        adj_[v].push_back(w);
    }

    /**
     * Makes a DFS Topological Sort of the complete graph and returns result vector of sorting algorithm.
     */
    vector<unsigned int> topological_sort();
};


// A recursive function used by topological_sort
void DfsTopoSort::topological_sort_util(unsigned int v, bool visited[],
                                stack<unsigned int>& stc)
{
    // Mark the current node as visited.
    visited[v] = true;

    // Recur for all the vertices
    // adjacent to this vertex
    list<unsigned int>::iterator i;
    for (i = adj_[v].begin(); i != adj_[v].end(); ++i)
        if (!visited[*i])
        	topological_sort_util(*i, visited, stc);

    // Push current vertex to stack
    // which stores result
    stc.push(v);
}

// The function to do Topological Sort.
// It uses recursive topological_sort_util()
vector<unsigned int> DfsTopoSort::topological_sort()
{
    stack<unsigned int> stc;

    // Mark all the vertices as not visited
    bool* visited = new bool[nv_];
    for (unsigned int i = 0; i < nv_; i++)
        visited[i] = false;

    // Call the recursive helper function
    // to store Topological
    // Sort starting from all
    // vertices one by one
    for (unsigned int i = 0; i < nv_; i++)
        if (visited[i] == false)
        	topological_sort_util(i, visited, stc);

    // Move contents of stack to vector
    vector<unsigned int> vec;
    vec.resize(stc.size());
    while (!stc.empty()) {
        vec[stc.size()-1] = stc.top();
        stc.pop();
    }

    return vec;
}

#endif /* DFS_TOPO_SORT_HH_ */
