/*
 * field_elementwise_impl.hh
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

/// Implementation.

namespace IT = Input::Type;

template <int spacedim, class Value>
it::Record FieldElementwise<spacedim, Value>::input_type
    = FieldElementwise<spacedim, Value>::get_input_type(FieldBase<spacedim, Value>::input_type, NULL);


template <int spacedim, class Value>
Input::Type::Record FieldElementwise<spacedim, Value>::get_input_type(
        Input::Type::AbstractRecord &a_type, typename Value::ElementInputType *eit
        )
{
    it::Record type=
        it::Record("FieldElementwise", FieldBase<spacedim,Value>::template_name()+" Field constant in space.")
        .derive_from(a_type)
        .declare_key("gmsh_file", IT::FileName::input(), IT::Default::obligatory(),
                "Input file with ASCII GMSH file format.")
        .declare_key("field_name", IT::String(), IT::Default::obligatory(),
                "The values of the Field are read from the $ElementData section with field name given by this key.")
        .close();

    return type;
}



template <int spacedim, class Value>
FieldElementwise<spacedim, Value>::FieldElementwise( unsigned int n_comp)
: FieldBase<spacedim, Value>(n_comp), data_(NULL), mesh_(NULL), reader_(NULL)
{
    n_components_ = this->value_.n_rows() * this->value_.n_cols();
}



template <int spacedim, class Value>
void FieldElementwise<spacedim, Value>::init_from_input(const Input::Record &rec) {
    FilePath input_file = rec.val<FilePath>("gmsh_file");
    ASSERT( reader_ == NULL, "Multiple call of init_from_input.\n");
    reader_ = new GmshMeshReader(input_file);

    field_name_ = rec.val<std::string>("field_name");
}



template <int spacedim, class Value>
void FieldElementwise<spacedim, Value>::set_data_row(unsigned int boundary_idx, typename Value::return_type &value) {
    Value ref(value);
    ASSERT( this->value_.n_cols() == ref.n_cols(), "Size of variable vectors do not match.\n" );
    typename Value::element_type *ptr=(typename Value::element_type *) ( data_+(bulk_size_+boundary_idx)*n_components_);
    for(unsigned int row=0; row < ref.n_rows(); row++)
        for(unsigned int col=0; col < ref.n_cols(); col++, ptr++)
            *ptr = ref(row,col);
}


template <int spacedim, class Value>
void FieldElementwise<spacedim, Value>::set_time(double time) {
    ASSERT(mesh_, "Null mesh pointer of elementwise field: %s, did you call set_mesh()?\n", field_name_.c_str());
    ASSERT(data_, "Null data pointer.\n");
    if (reader_ == NULL) return;

    GMSH_DataHeader search_header;
    search_header.actual=false;
    search_header.field_name=field_name_;
    search_header.n_components=n_components_;
    search_header.n_entities=n_entities_;
    search_header.time=time;

    reader_->read_element_data(search_header, data_, mesh_->all_elements_id() );
}



template <int spacedim, class Value>
void FieldElementwise<spacedim, Value>::set_mesh(Mesh *mesh) {
    // set mesh only once
    if (mesh_ != NULL) return;

    mesh_=mesh;
    bulk_size_=mesh_->element.size();
    n_entities_=mesh_->element.size() + mesh_->bc_elements.size();

    // allocate
    data_size_ = n_entities_ * n_components_;
    data_ = new double[data_size_];
    std::fill(data_, data_ + data_size_, 0.0);

}



/**
 * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
 */
template <int spacedim, class Value>
typename Value::return_type const & FieldElementwise<spacedim, Value>::value(const Point<spacedim> &p, const ElementAccessor<spacedim> &elm)
{

        unsigned int idx = elm.idx();
        if (elm.is_boundary()) idx +=bulk_size_;
        idx*=n_components_;

        return Value::from_raw(this->r_value_, (typename Value::element_type *)(data_+idx));
}



/**
 * Returns std::vector of scalar values in several points at once.
 */
template <int spacedim, class Value>
void FieldElementwise<spacedim, Value>::value_list (const std::vector< Point<spacedim> >  &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list)
{

    ASSERT_SIZES( point_list.size(), value_list.size() );
    if (boost::is_floating_point< typename Value::element_type>::value) {
        unsigned int idx = elm.idx();
        if (elm.is_boundary()) idx +=bulk_size_;
        idx*=n_components_;

        typename Value::return_type const &ref = Value::from_raw(this->r_value_, (typename Value::element_type *)(data_+idx));
        for(unsigned int i=0; i< value_list.size(); i++) value_list[i] = ref;
    } else {
        xprintf(UsrErr, "FieldElementwise is not implemented for discrete return types.\n");
    }
}



template <int spacedim, class Value>
FieldElementwise<spacedim, Value>::~FieldElementwise() {
    if (data_ != NULL) delete [] data_;
    if (reader_ != NULL) delete reader_;
}



#endif /* FIELD_ELEMENTWISE_IMPL_HH_ */
