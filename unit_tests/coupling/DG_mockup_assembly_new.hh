#ifndef DG_MOCKUP_ASSEMBLYA_HH_
#define DG_MOCKUP_ASSEMBLYA_HH_

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "DG_mockup.hh"
#include "fem/fe_p.hh"
#include "fem/patch_fe_values.hh"
#include "quadrature/quadrature_lib.hh"
#include "coupling/balance.hh"
#include "fem/element_cache_map.hh"
#include "transport/assembly_dg.hh"



template <unsigned int dim, class TEqFields, class TEqData>
class MassEvalFields : public MassAssemblyDG<dim, TEqFields, TEqData>
{
public:
    typedef TEqFields EqFields;
    typedef TEqData EqData;

    static constexpr const char * name() { return "MassAssemblyDG"; }

    /// Constructor.
    MassEvalFields(EqFields *eq_fields, EqData *eq_data, AssemblyInternals *asm_internals)
    : MassAssemblyDG<dim, TEqFields, TEqData>(eq_fields, eq_data, asm_internals) {}

    /// Destructor.
    ~MassEvalFields() {}

    void cell_integral(FMT_UNUSED DHCellAccessor cell, FMT_UNUSED unsigned int element_patch_idx) override {}

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


template <unsigned int dim, class TEqFields, class TEqData>
class StiffnessEvalFields : public StiffnessAssemblyDG<dim, TEqFields, TEqData>
{
public:
    typedef TEqFields EqFields;
    typedef TEqData EqData;

    static constexpr const char * name() { return "StiffnessAssemblyDG"; }

    /// Constructor.
    StiffnessEvalFields(EqFields *eq_fields, EqData *eq_data, AssemblyInternals *asm_internals)
    : StiffnessAssemblyDG<dim, TEqFields, TEqData>(eq_fields, eq_data, asm_internals) {}

    /// Destructor.
    ~StiffnessEvalFields() {}

    /// Assembles the cell (volume) integral into the stiffness matrix.
    inline void cell_integral(FMT_UNUSED DHCellAccessor cell, FMT_UNUSED unsigned int element_patch_idx) override {}

    /// Assembles the fluxes on the boundary.
    inline void boundary_side_integral(FMT_UNUSED DHCellSide cell_side) override {}


    /// Assembles the fluxes between sides of elements of the same dimension.
    inline void edge_integral(FMT_UNUSED RangeConvert<DHEdgeSide, DHCellSide> edge_side_range) override {}


    /// Assembles the fluxes between elements of different dimensions.
    inline void dimjoin_intergral(FMT_UNUSED DHCellAccessor cell_lower_dim, FMT_UNUSED DHCellSide neighb_side) override {}

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


template <unsigned int dim, class TEqFields, class TEqData>
class SourcesEvalFields : public SourcesAssemblyDG<dim, TEqFields, TEqData>
{
public:
    typedef TEqFields EqFields;
    typedef TEqData EqData;

    static constexpr const char * name() { return "SourcesAssemblyDG"; }

    /// Constructor.
    SourcesEvalFields(EqFields *eq_fields, EqData *eq_data, AssemblyInternals *asm_internals)
    : SourcesAssemblyDG<dim, TEqFields, TEqData>(eq_fields, eq_data, asm_internals) {}

    /// Destructor.
    ~SourcesEvalFields() {}

    void cell_integral(FMT_UNUSED DHCellAccessor cell, FMT_UNUSED unsigned int element_patch_idx) override {}

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


#endif /* DG_MOCKUP_ASSEMBLYA_HH_ */

