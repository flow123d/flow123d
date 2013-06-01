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
 * - move raw access resolution functions from FieldValues_ into FieldElementwise
 * - allow elementwise int or FieldEnum data with optimal storage buffer, this needs
 *   templated GMSH reader
 * - allow initialization of multiple fields by one reader
 * - allow common storage for more elementwise fields to have values for one element on one place
 */

#include "system/system.hh"
#include "fields/field_base.hh"

class GmshMeshReader;

template <int spacedim, class Value>
class FieldElementwise : public FieldBase<spacedim, Value>
{
public:

    FieldElementwise(unsigned int n_comp=0);

    /**
     * Temporary solution (as well as the whole this class before we use DofHandlers) how to
     * build FiledElementwise on an existing array of values on elements.
     */
    FieldElementwise(double *data_ptr, unsigned int n_components, unsigned int size );

    static Input::Type::Record input_type;

    static Input::Type::Record get_input_type(Input::Type::AbstractRecord &a_type, typename Value::ElementInputType *eit);

    virtual void init_from_input(const Input::Record &rec);

    /**
     * Set row of boundary data. Used to implement old BC input.
     */
    void set_data_row(unsigned int boundary_idx, typename Value::return_type &value);

    /**
     * Update time and possibly update data from GMSH file.
     */
    virtual bool set_time(double time);

    /**
     * Has to be set before calling init_from_input. This also
     * allocate and initialize internal buffer. Do nothing if mesh is already set.
     *
     * See also description of the FieldBase<...>::set_mesh.
     */
    virtual void set_mesh(Mesh *mesh, bool boundary_domain);


    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    virtual typename Value::return_type const &value(const Point<spacedim> &p, const ElementAccessor<spacedim> &elm);

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void value_list (const std::vector< Point<spacedim> >  &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list);


    virtual ~FieldElementwise();

private:
    /// If the data vector is provided at construction, we disallow initialization form input.
    bool allow_init_from_input;
    /**
     * Is set in set_mesh method. Value true means, that we accept only boundary element accessors in the @p value method.
     * TODO: temporary solution until we have separate mesh for the boundary part
     */
    bool boundary_domain_;
    /// number of bulk data lines in the buffer
    //unsigned int bulk_size_;
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
