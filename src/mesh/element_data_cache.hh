/*
 * element_data_cache.hh
 *
 *  Created on: Jan 28, 2013
 *      Author: jb
 */

#ifndef ELEMENT_DATA_CACHE_HH_
#define ELEMENT_DATA_CACHE_HH_

#include <mesh/mesh.h>


class ElementDataCache {
public:
	class QuantityStorageBase {

	};

	template <class T>
	class QuantityStorage : public QuantityStorageBase {
	public:
		typedef std::shared_ptr< std::vector<T> > ComponentData;
	    typedef unsigned int VectorSize; // number of T type values per element

	    QuantityStorage(Mesh &mesh, VectorSize vals_per_element)
	    : data_( std::make_shared< std::vector<T> >() ),
	      vals_per_element_(vals_per_element)
	    {
	    	data_.reserve(mesh.n_elements() * vals_per_element_);
	    }

	    ComponentData data_;
	    VectorSize vals_per_element_;
	};

    typedef std::map<std::string, QuantityStorageBase*> TimeStep;

};

#endif /* ELEMENT_DATA_CACHE_HH_ */
