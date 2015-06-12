/*
 * type_generic.hh
 *
 *  Created on: May 1, 2012
 *      Author: jb
 */

#ifndef TYPE_GENERIC_HH_
#define TYPE_GENERIC_HH_



#include <input/type_base.hh>



namespace Input {

namespace Type {


/**
 * Class for representing Record key of any type.
 *
 * Instances of this class are used only in generic types and during generation
 * of Record are replaced by types of IST (Integer, String, Selection etc.)
 */
class Parameter : public TypeBase {
public:
	Parameter(const string & parameter_name);

protected:
	const string name_;
};



} // closing namespace Type
} // closing namespace Input



#endif /* TYPE_GENERIC_HH_ */
