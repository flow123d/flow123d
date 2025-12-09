#ifndef ELASTICITY_MOCKUP_ASSEMBLYA_HH_
#define ELASTICITY_MOCKUP_ASSEMBLYA_HH_

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "mechanics/assembly_elasticity.hh"
#include "elasticity_mockup.hh"
#include "fem/patch_fe_values.hh"


template <unsigned int dim, class TEqData>
class StiffnessEvalFields : public StiffnessAssemblyElasticity<dim, TEqData>
{
public:
    typedef typename TEqData::EqFields EqFields;
    typedef TEqData EqData;

    static constexpr const char * name() { return "StiffnessAssemblyElasticity"; }

    /// Constructor.
    StiffnessEvalFields(EqData *eq_data, PatchFEValues<3> *fe_values)
    : StiffnessAssemblyElasticity<dim, TEqData>(eq_data, fe_values) {}

    /// Destructor.
    ~StiffnessEvalFields() {}

    /// Assembles the cell (volume) integral into the stiffness matrix.
    inline void cell_integral(FMT_UNUSED DHCellAccessor cell, FMT_UNUSED unsigned int element_patch_idx) override {}

    /// Assembles the fluxes on the boundary.
    inline void boundary_side_integral(FMT_UNUSED DHCellSide cell_side) override {}

    /// Assembles the fluxes between elements of different dimensions.
    inline void dimjoin_intergral(FMT_UNUSED DHCellAccessor cell_lower_dim, FMT_UNUSED DHCellSide neighb_side) override {}

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


template <unsigned int dim, class TEqData>
class RhsEvalFields : public RhsAssemblyElasticity<dim, TEqData>
{
public:
    typedef typename TEqData::EqFields EqFields;
    typedef TEqData EqData;

    static constexpr const char * name() { return "RhsAssemblyElasticity"; }

    /// Constructor.
    RhsEvalFields(EqData *eq_data, PatchFEValues<3> *fe_values)
    : RhsAssemblyElasticity<dim, TEqData>(eq_data, fe_values) {}

    /// Destructor.
    ~RhsEvalFields() {}

    /// Assembles the cell (volume) integral into the stiffness matrix.
    inline void cell_integral(FMT_UNUSED DHCellAccessor cell, FMT_UNUSED unsigned int element_patch_idx) override {}

    /// Assembles the fluxes on the boundary.
    inline void boundary_side_integral(FMT_UNUSED DHCellSide cell_side) override {}

    /// Assembles the fluxes between elements of different dimensions.
    inline void dimjoin_intergral(FMT_UNUSED DHCellAccessor cell_lower_dim, FMT_UNUSED DHCellSide neighb_side) override {}

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


#endif /* ELASTICITY_MOCKUP_ASSEMBLYA_HH_ */

