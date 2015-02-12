/*
 * field_elementwise.impl.hh
 *
 *  Created on: Jan 23, 2013
 *      Author: jb
 */

#ifndef FIELD_ELEMENTWISE_IMPL_HH_
#define FIELD_ELEMENTWISE_IMPL_HH_


#include "fields/field_elementwise.hh"
#include "system/file_path.hh"
#include "input/input_type.hh"
#include "mesh/msh_gmshreader.h"
#include "mesh/reader_instances.hh"

/// Implementation.

namespace IT = Input::Type;

template <int spacedim, class Value>
it::Record FieldElementwise<spacedim, Value>::input_type
    = FieldElementwise<spacedim, Value>::get_input_type(FieldAlgorithmBase<spacedim, Value>::input_type, NULL);


template <int spacedim, class Value>
Input::Type::Record FieldElementwise<spacedim, Value>::get_input_type(
        Input::Type::AbstractRecord &a_type, const typename Value::ElementInputType *eit
        )
{
    it::Record type=
        it::Record("FieldElementwise", FieldAlgorithmBase<spacedim,Value>::template_name()+" Field constant in space.")
        .derive_from(a_type)
        .declare_key("gmsh_file", IT::FileName::input(), IT::Default::obligatory(),
                "Input file with ASCII GMSH file format.")
        .declare_key("field_name", IT::String(), IT::Default::obligatory(),
                "The values of the Field are read from the \\$ElementData section with field name given by this key.")
        .close();

    return type;
}


template <int spacedim, class Value>
const int FieldElementwise<spacedim, Value>::registrar =
		Input::register_class< FieldElementwise<spacedim, Value>, unsigned int >("FieldElementwise");



template <int spacedim, class Value>
FieldElementwise<spacedim, Value>::FieldElementwise( unsigned int n_comp)
: FieldAlgorithmBase<spacedim, Value>(n_comp),
  internal_raw_data(true), data_(NULL), mesh_(NULL)

{
    n_components_ = this->value_.n_rows() * this->value_.n_cols();
}



template <int spacedim, class Value>
FieldElementwise<spacedim, Value>::FieldElementwise(double *data_ptr, unsigned int n_components, unsigned int size )
: FieldAlgorithmBase<spacedim, Value>(n_components),
  internal_raw_data(false), data_size_(size), data_(data_ptr), mesh_(NULL)
{
    n_components_ = this->value_.n_rows() * this->value_.n_cols();
}



template <int spacedim, class Value>
void FieldElementwise<spacedim, Value>::init_from_input(const Input::Record &rec) {
	cout << string(reader_file_) << endl;
    ASSERT( internal_raw_data, "Trying to initialize internal FieldElementwise from input.");
    ASSERT( reader_file_ == FilePath(), "Multiple call of init_from_input.\n");
    reader_file_ = FilePath( rec.val<FilePath>("gmsh_file") );
    ReaderInstances::instance()->get_reader(reader_file_);

    field_name_ = rec.val<std::string>("field_name");
}



template <int spacedim, class Value>
void FieldElementwise<spacedim, Value>::set_data_row(unsigned int boundary_idx, typename Value::return_type &value) {
    Value ref(value);
    ASSERT( this->value_.n_cols() == ref.n_cols(), "Size of variable vectors do not match.\n" );
    ASSERT( mesh_, "Null mesh pointer of elementwise field: %s, did you call set_mesh()?\n", field_name_.c_str());
    ASSERT( boundary_domain_ , "Method set_data_row can be used only for boundary fields.");
    typename Value::element_type *ptr=(typename Value::element_type *) ( data_+(boundary_idx)*n_components_);
    for(unsigned int row=0; row < ref.n_rows(); row++)
        for(unsigned int col=0; col < ref.n_cols(); col++, ptr++)
            *ptr = ref(row,col);
}


template <int spacedim, class Value>
bool FieldElementwise<spacedim, Value>::set_time(double time) {
    ASSERT(mesh_, "Null mesh pointer of elementwise field: %s, did you call set_mesh()?\n", field_name_.c_str());
    ASSERT(data_, "Null data pointer.\n");
    if ( reader_file_ == FilePath() ) return false;

    //walkaround for the steady time governor - there is no data to be read in time==infinity
    //TODO: is it possible to check this before calling set_time?
    if (time == numeric_limits< double >::infinity()) return false;
    
    GMSH_DataHeader search_header;
    search_header.actual=false;
    search_header.field_name=field_name_;
    search_header.n_components=n_components_;
    search_header.n_entities=n_entities_;
    search_header.time=time;


    ReaderInstances::instance()->get_reader(reader_file_)->read_element_data(search_header, data_, mesh_->elements_id_maps(boundary_domain_) );
    return search_header.actual;
}



template <int spacedim, class Value>
void FieldElementwise<spacedim, Value>::set_mesh(const Mesh *mesh, bool boundary_domain) {
    // set mesh only once or to same value
    ASSERT(mesh_ == nullptr || mesh_ == mesh, "Trying to change mesh of the FieldElementwise.");
    boundary_domain_ = boundary_domain;

    mesh_=mesh;
    if (boundary_domain_) {
        n_entities_=mesh_->bc_elements.size();
    } else {
        n_entities_=mesh_->n_elements();
    }

    // allocate
    if (data_ == NULL) {
        data_size_ = n_entities_ * n_components_;
        data_ = new double[data_size_];
        std::fill(data_, data_ + data_size_, 0.0);
    }

}



/**
 * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
 */
template <int spacedim, class Value>
typename Value::return_type const & FieldElementwise<spacedim, Value>::value(const Point &p, const ElementAccessor<spacedim> &elm)
{
        ASSERT( elm.is_elemental(), "FieldElementwise works only for 'elemental' ElementAccessors.\n");
        ASSERT( elm.is_boundary() == boundary_domain_, "Trying to get value of FieldElementwise '%s' for wrong ElementAccessor type (boundary/bulk).\n", field_name_.c_str() );

        unsigned int idx = n_components_*elm.idx();

        return Value::from_raw(this->r_value_, (typename Value::element_type *)(data_+idx));
}



/**
 * Returns std::vector of scalar values in several points at once.
 */
template <int spacedim, class Value>
void FieldElementwise<spacedim, Value>::value_list (const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list)
{
    ASSERT( elm.is_elemental(), "FieldElementwise works only for 'elemental' ElementAccessors.\n");
    ASSERT( elm.is_boundary() == boundary_domain_, "Trying to get value of FieldElementwise '%s' for wrong ElementAccessor type (boundary/bulk).\n", field_name_.c_str() );
    ASSERT_EQUAL( point_list.size(), value_list.size() );
    if (boost::is_floating_point< typename Value::element_type>::value) {
        unsigned int idx = n_components_*elm.idx();

        typename Value::return_type const &ref = Value::from_raw(this->r_value_, (typename Value::element_type *)(data_+idx));
        for(unsigned int i=0; i< value_list.size(); i++) {
            ASSERT( Value(value_list[i]).n_rows()==this->value_.n_rows(),
                    "value_list[%d] has wrong number of rows: %d; should match number of components: %d\n",
                    i, Value(value_list[i]).n_rows(),this->value_.n_rows());

            value_list[i] = ref;
        }
    } else {
        xprintf(UsrErr, "FieldElementwise is not implemented for discrete return types.\n");
    }
}



template <int spacedim, class Value>
FieldElementwise<spacedim, Value>::~FieldElementwise() {
    if (internal_raw_data && data_ != NULL) delete [] data_;
}



#endif /* FIELD_ELEMENTWISE_IMPL_HH_ */
