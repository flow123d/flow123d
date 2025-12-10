#ifndef DG_MOCKUP_ASSEMBLYA_HH_
#define DG_MOCKUP_ASSEMBLYA_HH_

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "DG_mockup.hh"
#include "fem/fe_p.hh"
#include "fem/patch_fe_values.hh"
#include "fem/op_factory.hh"
#include "quadrature/quadrature_lib.hh"
#include "coupling/balance.hh"
#include "fem/element_cache_map.hh"
#include "transport/assembly_dg.hh"



template <unsigned int dim, class TEqData>
class MassEvalFields : public MassAssemblyDG<dim, TEqData>
{
public:
    typedef typename TEqData::EqFields EqFields;
    typedef TEqData EqData;

    static constexpr const char * name() { return "DG_Mass_Assembly"; }

    /// Constructor.
    MassEvalFields(EqData *eq_data, PatchFEValues<3> *fe_values)
    : MassAssemblyDG<dim, TEqData>(eq_data, fe_values) {}

    /// Destructor.
    ~MassEvalFields() {}

    void cell_integral(FMT_UNUSED DHCellAccessor cell, FMT_UNUSED unsigned int element_patch_idx) override {}

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


template <unsigned int dim, class TEqData>
class StiffnessEvalFields : public StiffnessAssemblyDG<dim, TEqData>
{
public:
    typedef typename TEqData::EqFields EqFields;
    typedef TEqData EqData;

    static constexpr const char * name() { return "DG_Stiffness_Assembly"; }

    /// Constructor.
    StiffnessEvalFields(EqData *eq_data, PatchFEValues<3> *fe_values)
    : StiffnessAssemblyDG<dim, TEqData>(eq_data, fe_values) {}

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


template <unsigned int dim, class TEqData>
class SourcesEvalFields : public SourcesAssemblyDG<dim, TEqData>
{
public:
    typedef typename TEqData::EqFields EqFields;
    typedef TEqData EqData;

    static constexpr const char * name() { return "DG_Sources_Assembly"; }

    /// Constructor.
    SourcesEvalFields(EqData *eq_data, PatchFEValues<3> *fe_values)
    : SourcesAssemblyDG<dim, TEqData>(eq_data, fe_values) {}

    /// Destructor.
    ~SourcesEvalFields() {}

    void cell_integral(FMT_UNUSED DHCellAccessor cell, FMT_UNUSED unsigned int element_patch_idx) override {}

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


#endif /* DG_MOCKUP_ASSEMBLYA_HH_ */

