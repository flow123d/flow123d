#ifndef DARCY_MOCKUP_ASSEMBLYA_HH_
#define DARCY_MOCKUP_ASSEMBLYA_HH_

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "darcy_mockup.hh"
#include "fem/fe_p.hh"
#include "fem/patch_fe_values.hh"
#include "fem/op_factory.hh"
#include "quadrature/quadrature_lib.hh"
#include "coupling/balance.hh"
#include "fem/element_cache_map.hh"
#include "flow/assembly_lmh.hh"



template <unsigned int dim, class TEqFields, class TEqData>
class MHMatrixEvalFields : public MHMatrixAssemblyLMH<dim, TEqFields, TEqData>
{
public:
    typedef TEqFields EqFields;
    typedef TEqData EqData;

    static constexpr const char * name() { return "MHMatrixAssemblyLMH"; }

    /// Constructor.
    MHMatrixEvalFields(EqFields *eq_fields, EqData *eq_data)
    : MHMatrixAssemblyLMH<dim, TEqFields, TEqData>(eq_fields, eq_data) {}

    /// Destructor.
    ~MHMatrixEvalFields() {}

    /// Assembles the cell (volume) integral into the stiffness matrix.
    inline void cell_integral(FMT_UNUSED DHCellAccessor cell, FMT_UNUSED unsigned int element_patch_idx) override {}

    /// Assembles the fluxes on the boundary.
    inline void boundary_side_integral(FMT_UNUSED DHCellSide cell_side) override {}

    /// Assembles the fluxes between elements of different dimensions.
    inline void dimjoin_intergral(FMT_UNUSED DHCellAccessor cell_lower_dim, FMT_UNUSED DHCellSide neighb_side) override {}

    /// Implements @p AssemblyBase::begin.
    void begin() override {}

    /// Implements @p AssemblyBase::end.
    void end() override {}

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


#endif /* DARCY_MOCKUP_ASSEMBLYA_HH_ */

