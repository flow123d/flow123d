/*
 * multi_field.impl.hh
 *
 *  Created on: Feb 13, 2014
 *      Author: jb
 */

#ifndef MULTI_FIELD_IMPL_HH_
#define MULTI_FIELD_IMPL_HH_


#include "multi_field.hh"
#include "fields/field_algo_base.hh"

namespace it = Input::Type;

/******************************************************************************************
 * Implementation of MultiField<...>
 */

template<int spacedim, class Value>
MultiField<spacedim, Value>::MultiField()
: FieldCommon()
{}



template<int spacedim, class Value>
void MultiField<spacedim, Value>::init( const vector<string> &names) {
    sub_fields_.resize( names.size() );
    sub_names_ = names;
    for(unsigned int i_comp=0; i_comp < size(); i_comp++)
    {
    	sub_fields_[i_comp].units( units() );

    	if (sub_names_[i_comp].length() == 0)
    		sub_fields_[i_comp].name( name() );
    	else
    		sub_fields_[i_comp].name( sub_names_[i_comp] + "_" + name());
    }
}



template<int spacedim, class Value>
it::Record &  MultiField<spacedim,Value>::get_input_type() {
	static it::Record type= it::Record("MultiField", "Record for all time-space functions.")
		.has_obligatory_type_key()
		.declare_key("component_names", it::Array( it::String() ), it::Default::read_time("Can be get from source of MultiField."),
			"Names of MultiField components.")
		.declare_key("common", transposed_field_.get_input_type(), it::Default::optional(),
			"Supplied to components subtree.")
		.declare_key("components", it::Array( sub_field_type_.get_input_type() ), it::Default::read_time("Converts from 'common' key."),
			"Components of Multifield.");

	return type;
}


template<int spacedim, class Value>
void MultiField<spacedim, Value>::set_limit_side(LimitSide side)
{
	for ( SubFieldType &field : sub_fields_)
		field.set_limit_side(side);
}


template<int spacedim, class Value>
bool MultiField<spacedim, Value>::set_time(
		const TimeGovernor &time)
{
	bool any=false;
	for( SubFieldType &field : sub_fields_) {
		if (field.set_time(time))
			any = true;
	}
    return any;
}



template<int spacedim, class Value>
void MultiField<spacedim, Value>::set_mesh(const Mesh &mesh) {
    shared_->mesh_ = &mesh;
    for(unsigned int i_comp=0; i_comp < size(); i_comp++)
        sub_fields_[i_comp].set_mesh(mesh);
}


template<int spacedim, class Value>
void MultiField<spacedim, Value>::copy_from(const FieldCommon & other) {
	if (typeid(other) == typeid(*this)) {
		auto  const &other_field = dynamic_cast<  MultiField<spacedim, Value> const &>(other);
		this->operator=(other_field);
	} else if (typeid(other) == typeid(SubFieldType)) {
		auto  const &other_field = dynamic_cast<  SubFieldType const &>(other);
		sub_fields_.resize(1);
		sub_fields_[0] = other_field;
	}
}



template<int spacedim, class Value>
void MultiField<spacedim, Value>::output(OutputTime *stream)
{
	// currently we cannot output boundary fields
	if (!is_bc())
		stream->register_data(this->output_type(), *this);
}





template<int spacedim, class Value>
bool MultiField<spacedim, Value>::is_constant(Region reg) {
	bool const_all=false;
	for(auto field : sub_fields_) const_all = const_all || field.is_constant(reg);
	return const_all;
}



template<int spacedim, class Value>
typename Field<spacedim,Value>::FieldBasePtr MultiField<spacedim, Value>::MultiFieldFactory::create_field(Input::Record rec, const FieldCommon &field) {
	Input::Iterator<Input::AbstractRecord> it_common = rec.find<Input::AbstractRecord>("common");
	Input::Iterator<Input::Array> it_components = rec.find<Input::Array>("components");
	if (it_common && !it_components) {
		it_common->transpose_to( rec, "components", spacedim );
		it_components = rec.find<Input::Array>("components");
		for (auto it = it_components->begin<Input::AbstractRecord>(); it != it_components->end(); ++it)
		{
			typename Field<spacedim,Value>::FieldBasePtr func =
					Field<spacedim,Value>::FieldBaseType::function_factory( (*it), field.n_comp() );
			//TODO: set func to ??
		}
	}
	return NULL;
}



#endif /* MULTI_FIELD_IMPL_HH_ */
