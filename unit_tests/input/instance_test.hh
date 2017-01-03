/*
 * instance_test.h
 *
 *  Created on: Apr 4, 2012
 *      Author: jb
 *

 *
 */

#ifndef INSTANCE_TEST_HH_
#define INSTANCE_TEST_HH_


namespace it = Input::Type;

#include <input/input_type.hh>

std::string get_tmp_name() {
	static std::vector<std::string> names = { "Integer", "Double", "String" };
	static int idx = 0;
	return names[idx++];

}

template<class Type>
class InstanceTest
{
public:
    /**
     * Returns input type Record.
     */
    static Input::Type::Record & get_input_type(std::string tmp_name)
    {
        return it::Record("GenericRecord", "Test generic record of "+tmp_name)
            .root_of_generic_subtree()
            .declare_key("a_val", Parameter("a_val"), Default::obligatory(), "value A" )
            .declare_key("b_val", Parameter("b_val"), Default::obligatory(), "value B" )
            .declare_key("description", it::String(), Default::obligatory(), "description" )
    		.close();

    }

    /**
     * Returns parameterized input type Record.
     */
    static const Input::Type::Instance & get_input_type_instance(std::string tmp_name)
    {
    	std::vector<it::TypeBase::ParameterPair> param_vec;
    	param_vec.push_back( std::make_pair("a_val", std::make_shared<Type>()) );
    	param_vec.push_back( std::make_pair("b_val", std::make_shared<Type>()) );

    	return it::Instance(get_input_type(tmp_name), param_vec).close();
    }
private:
    static const int registrar;
};

template <class Type>
const int InstanceTest<Type>::registrar = InstanceTest<Type>::get_input_type( get_tmp_name() ).size();


#endif /* INSTANCE_TEST_HH_ */
