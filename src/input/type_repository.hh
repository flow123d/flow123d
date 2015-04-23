/*
 * type_repository.hh
 *
 *  Created on: May 7, 2012
 *      Author: jb
 */

#ifndef TYPE_REPOSITORY_HH_
#define TYPE_REPOSITORY_HH_


#include <map>
#include "type_base.hh"


namespace Input {


/**
 * The Singleton class TypeRepository serves for handling the lazy-evaluated input types, derived from the base class
 * LazyType. When all static variables are initialized, the method TypeRepository::instance().finish() can be called
 * in order to finish initialization of lazy types such as Records, AbstractRecords, Arrays and Selections.
 * Selections have to be finished after all other types since they are used by AbstractRecords to register all
 * derived types. For this reason TypeRepository contains two arrays - one for Selections, one for the rest.
 *
 * This is list of unique instances that may contain raw pointers to possibly not yet constructed
 * (static) objects. Unique instance is the instance that creates unique instance of the data class in pimpl idiom.
 * These has to be completed/finished before use.
 *
 */
template <class T>
class TypeRepository {
public:
	static_assert(std::is_base_of<Type::TypeBase, T>::value,
	        "T must be a descendant of Input::Type::TypeBase"
	    );

	/**
     * Type of hash values used in associative array that translates key names to indices in Record and Selection.
     *
     * For simplicity, we currently use whole strings as "hash".
     */
    typedef string KeyHash;

    typedef std::map< KeyHash, boost::shared_ptr<T> > TypeRepositoryMap;

    /// Hash function.
    inline static KeyHash key_hash(const string &str) {
        return (str);
    }

    static boost::shared_ptr<T> add_type(const T & type);

};


template <class T>
boost::shared_ptr<T> TypeRepository<T>::add_type(const T & type) {
    static TypeRepositoryMap type_repository_map;
	KeyHash hash = TypeRepository<T>::key_hash(type.full_type_name());

	auto search = type_repository_map.find(hash);
	if (search != type_repository_map.end()) {
		return search->second;
	} else {
		auto type_ptr = boost::make_shared<T>( type );
		type_repository_map.insert( std::pair<KeyHash, boost::shared_ptr<T>>(hash,type_ptr) );
		return type_ptr;
	}
}

} // namespace Input


#endif /* TYPE_REPOSITORY_HH_ */
