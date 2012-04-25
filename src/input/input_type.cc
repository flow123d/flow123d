/*
 * input_type.cc
 *
 *  Created on: Mar 29, 2012
 *      Author: jb
 */


#include "input_type.hh"

namespace Input {
namespace Type {

std::ostream& operator<<(std::ostream& stream, const TypeBase& type) {
    return type.documentation(stream);
}

} // closing namespace Type
} // closing namespace Input


