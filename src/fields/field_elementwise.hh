/*
 * field_elementwise.hh
 *
 *  Created on: Jan 22, 2013
 *      Author: jb
 */

#ifndef FIELD_ELEMENTWISE_HH_
#define FIELD_ELEMENTWISE_HH_

/**
 * The simplest approximative Field. In future we want to replace it by pair: DofHandler + vector of dof values.
 *
 * In current setting:
 * This should implement mapping: (ElementAccessor.INDEX) -> value
 * Currently INDEX is global index of element (any dimension, both bulk and boundary part). In future INDEX should
 * probably consist of MESH_LEVEL (unique mesh, dimension, bulk/boundary, level of refinement) and INDEX within this level.
 * We want to make memory optimization since usually field lives either on boundary or on bulk part and some bulk fields live only on some dimension(s).
 * This can be achieved by two level indirection table of mesh_levelscontaining tables for indexes. We should test if the performance penalty is not to big.
 *
 * Currently, we just use one vector for bulk and one for boundary elements.
 *
 * TODO:
 * - input type
 * - init form input - read either from GMSH or from text file
 * - value
 *
 * - ve FieldValues implement inline constructor:
 *   FieldValues_<...>( return_type &, double * raw_data)
 *   for double set internal reference to *raw_data
 *   for arma types set reference to return_type and set its internal pointer mem to raw_data
 *
 *   carfully test !!
 *
 *   this allows as store raw data while providing fast access
 * - next solve initialization problem
 *   reader:
 *   while header.time < time .. read data, save header data in GMSH reader
 */

#include "system/system.hh"
#include "fields/field_base.hh"

class GmshMeshReader;

template <int spacedim, class Value>
class FieldElementwise : public FieldBase<spacedim, Value>
{
public:

    FieldElementwise(unsigned int n_comp=0);

    static Input::Type::Record input_type;

    static Input::Type::Record get_input_type(Input::Type::AbstractRecord &a_type, typename Value::ElementInputType *eit);

    virtual void init_from_input(const Input::Record &rec);

    /**
     * Update time and possibly update data from GMSH file.
     */
    virtual void set_time(double time);

    /**
     * Has to be set before calling init_from_input.
     */
    virtual void set_mesh(Mesh *mesh);

    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    virtual typename Value::return_type const &value(const Point<spacedim> &p, ElementAccessor<spacedim> &elm);

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void value_list (const std::vector< Point<spacedim> >  &point_list, ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list);


    virtual ~FieldElementwise();

private:
    /// number of bulk data lines in the buffer
    unsigned int bulk_size_;
    /// Allocated size of data_ buffer
    unsigned int data_size_;
    /// Raw buffer of n_entities rows each containing Value::size() doubles.
    double *data_;
    /// Number of rows in @p data_ buffer.
    unsigned int n_entities_;
    /// Size of Value
    unsigned int n_components_;

    GmshMeshReader *reader_;
    Mesh *mesh_;
    std::string field_name_;
};


#endif /* FIELD_ELEMENTWISE_HH_ */
