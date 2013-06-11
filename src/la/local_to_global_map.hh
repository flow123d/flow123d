/*
 * local_to_global_map.hh
 *
 *  Created on: Mar 9, 2012
 *      Author: jb
 */

#ifndef LOCAL_TO_GLOBAL_MAP_HH_
#define LOCAL_TO_GLOBAL_MAP_HH_

#include <set>
#include <vector>

/// Using Boost shared pointer.
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/smart_ptr/make_shared.hpp>

#include "system/system.hh"
#include "system/global_defs.h"

//namespace boost{    
//    template<class T> class shared_ptr;
//}

class Distribution;

/**
 * @brief class to manage local indices on sub-domain to global indices on domain
 *
 * Currently (March 2012) the local to global maps are managed by individual equations (e.g. DarcyFlow has el_4_loc .. there it is complicated
 * by different local indexing of elements, sides and edges). In fact el_4_loc etc. are new_local to old_global maps.
 * This map is created as follows:
 *
 * PETSC solver:
 * 1) create graph
 * 2) make partitioning
 * 2.5) create id_4_old (just used to create new_4_id, id was non continuous index not used anymore)
 * 3) call id_maps which:
 * 4) call ISPartitioninToNumbering : assign partitions to processors (identity, no optimization); make "distribution";
 *    make mapping (array of ints): old_local to new_global
 * 5) from this we crate AO (application ordering from PETSc)
 * 6) use AO to map identity array to old numbering -> creates map: new_global to old_global
 * 7) go through new_local continuous part and map it to old_global index
 * ---
 * So this produce new_local to old_global mapping and only for continuous part of local indices (not for ghost indices)
 * See that this works without any information about connectivity.
 *
 * METIS/BDDC solver:
 * 1) graph, partition of elements, maps for elements (no overlap, no ghost values) using id_maps function
 * 2) same for edges -> produce non overlapping local to global maps.
 * 3) use std::set to collect all dofs on local elements (pass through mesh) complexity: n*log(n) n-local number of dofs
 * 4) copy set to vector - makes local ordering arbitrary
 *
 * Future usage is:
 * - in mesh to map local idx of entities to global indices
 * - in dof handler to map local dofs idx to global
 *
 * In both cases the mapping is created by adding all global numbers on local subdomain. In order to keep local part of the map
 * continuous, we need Distribution that describes splitting of global indices into continuous blocks. This Distribution is known
 * as soon as we assign partitions to processors. So we assume that it is known at construction time.
 *
 */



class LocalToGlobalMap {
public:
   /**
    * Constructor. Starts filling of the map.
    *
    * @param distr Non overlapping distribution of global indices to processors in continuous blocks.
    * Local block of indices forms first part of the mapping and then nonlocal indices follows.
    * This constructor makes a deep copy of the distribution.
    *
    */
   LocalToGlobalMap(const Distribution &distr);

   /**
    * Same as the previous constructor, but just takes copy of shared pointer. This assumes that Distribution is allocated on the heap
    * by something like:
    *
    * boost::smart_ptr<Distribution> distr(new Distribution(...));
    */
   LocalToGlobalMap(boost::shared_ptr<Distribution> distr);

   /**
    * Insert a global index to the mapping.
    */
   void insert(const unsigned int idx);
   /**
    * Insert more indices at once.
    */
   void insert(const std::vector<unsigned int> &indices);
   /**
    * Finish filling stage. Creates mapping array, then you can use () operator to map indices.
    */
   void finalize();
   /**
    * Maps local index to the global one. For DEBUG, performs check for dimension.
    */
   inline unsigned int operator[] (const unsigned int local_idx) const
       {
           ASSERT_LESS( local_idx, global_indices_.size() );
           return global_indices_[local_idx];
       }

   /**
    * Returns size of local map.
    */
   unsigned int size() const
       { return global_indices_.size(); }

   /**
    * Returns smart_ptr to the Distribution of the global indices. Allow share this among many objects.
    */
   boost::shared_ptr<Distribution> &get_distr()
        { return distr_; }

   /**
    * Returns inner vector.
    */
   const std::vector<unsigned int> &get_map_vector() const
        { return global_indices_; }



private:
   /// auxiliary set of non-local part of the map
   std::set<unsigned int> *nonlocal_indices_;
   /// distribution of the global indices
   boost::shared_ptr<Distribution> distr_;
   /// mapping for all local indices
   std::vector<unsigned int> global_indices_;
};


#endif /* LOCAL_TO_GLOBAL_MAP_HH_ */
