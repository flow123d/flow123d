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
{
	this->multifield_ = true;
}



template<int spacedim, class Value>
const it::Instance &  MultiField<spacedim,Value>::get_input_type() {
	ASSERT(false, "This method can't be used for MultiField");

	static it::AbstractRecord abstract = it::AbstractRecord();
	static it::Instance inst = it::Instance( abstract, std::vector<it::TypeBase::ParameterPair>() );
	return inst;
}


template<int spacedim, class Value>
it::Record &  MultiField<spacedim,Value>::get_multifield_input_type() {
	static it::Record type= it::Record("MultiField", "Record for all time-space functions.")
		.has_obligatory_type_key()
		.declare_key("component_names", it::Array( it::String() ), it::Default::read_time("Can be get from source of MultiField."),
			"Names of MultiField components.")
		.declare_key("common", transposed_field_.get_input_type(), it::Default::optional(),
			"Supplied to components subtree.")
		.declare_key("components", it::Array( sub_field_type_.get_input_type() ), it::Default::read_time("Converts from 'common' key."),
			"Components of Multifield.")
		.close();

	return type;
}


template<int spacedim, class Value>
bool MultiField<spacedim, Value>::set_time(
		const TimeStep &time)
{
	// initialization of Multifield for first call
	if (sub_fields_.size() == 0) {
	    set_up_components();
	}

	// set time for sub fields
	bool any=false;
	for( SubFieldType &field : sub_fields_) {
		if (field.set_time(time))
			any = true;
	}
    return any;
}



template<int spacedim, class Value>
void MultiField<spacedim, Value>::set_mesh(const Mesh &mesh) {
	// test if mesh is not set
	if (shared_->mesh_ && shared_->mesh_ != &mesh) {
		THROW(ExcFieldMeshDifference() << EI_Field(name()) );
	}

    shared_->mesh_ = &mesh;
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
void MultiField<spacedim, Value>::output(std::shared_ptr<OutputTime> stream)
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
void MultiField<spacedim, Value>::set_up_components() {
	ASSERT(this->shared_->comp_names_.size(), "Vector of component names is empty!\n");
	ASSERT(this->shared_->mesh_, "Mesh is not set!\n");

    sub_fields_.resize( this->shared_->comp_names_.size() );
    for(unsigned int i_comp=0; i_comp < size(); i_comp++)
    {
    	sub_fields_[i_comp].units( units() );
    	sub_fields_[i_comp].set_mesh( *(shared_->mesh_) );
    	sub_fields_[i_comp].set_limit_side(this->limit_side_);
    	sub_fields_[i_comp].add_factory( std::make_shared<MultiFieldFactory>(i_comp) );
    	sub_fields_[i_comp].set_component_index(i_comp);

    	if (this->shared_->comp_names_[i_comp].length() == 0)
    		sub_fields_[i_comp].name( name() );
    	else {
    		sub_fields_[i_comp].name_ = this->shared_->comp_names_[i_comp] + "_" + name();
    		sub_fields_[i_comp].shared_->input_name_ = name();
    	}

    	sub_fields_[i_comp].set_input_list(this->shared_->input_list_);
    }
}



template<int spacedim, class Value>
typename Field<spacedim,Value>::FieldBasePtr MultiField<spacedim, Value>::MultiFieldFactory::create_field(Input::Record descriptor_rec, const FieldCommon &field) {
	Input::Record multifield_rec;
	if (descriptor_rec.opt_val(field.input_name(), multifield_rec));
	Input::Iterator<Input::AbstractRecord> it_common = multifield_rec.find<Input::AbstractRecord>("common");
	Input::Iterator<Input::Array> it_components = multifield_rec.find<Input::Array>("components");
	if (it_common && !it_components) {
		it_common->transpose_to( multifield_rec, "components", spacedim );
		it_components = multifield_rec.find<Input::Array>("components");
	}

	ASSERT(it_components, "Failed to fill 'components' array of multifield: %s.", field.input_name().c_str());
	ASSERT(index_ < it_components->size(), "Index of MultiField component is out of range.\n");

	unsigned int position = 0;
	for (auto it = it_components->begin<Input::AbstractRecord>(); it != it_components->end(); ++it, ++position)
	{
		if (index_ == position) {
			typename Field<spacedim,Value>::FieldBasePtr field_algo_base = Field<spacedim,Value>::FieldBaseType::function_factory( (*it), field.n_comp() );
			field_algo_base->set_component_idx(index_);
			return field_algo_base;
		}
	}
	return NULL;
}



#endif /* MULTI_FIELD_IMPL_HH_ */
