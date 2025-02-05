#ifndef ELASTICITY_MOCKUP_ASSEMBLYA_HH_
#define ELASTICITY_MOCKUP_ASSEMBLYA_HH_

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "elasticity_mockup.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "quadrature/quadrature_lib.hh"
#include "coupling/balance.hh"
#include "fields/field_value_cache.hh"


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim>
class Stiffness_FullAssembly : public AssemblyBase<dim>
{
public:
    typedef typename equation_data::EqFields EqFields;
    typedef typename equation_data::EqData EqData;

    static constexpr const char * name() { return "StiffnessAssembly"; }

    /// Constructor.
    Stiffness_FullAssembly(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(1), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = (ActiveIntegrals::bulk | ActiveIntegrals::coupling | ActiveIntegrals::boundary);
        this->used_fields_ += eq_fields_->cross_section;
        this->used_fields_ += eq_fields_->lame_mu;
        this->used_fields_ += eq_fields_->lame_lambda;
        this->used_fields_ += eq_fields_->dirichlet_penalty;
        this->used_fields_ += eq_fields_->bc_type;
        this->used_fields_ += eq_fields_->fracture_sigma;
    }

    /// Destructor.
    ~Stiffness_FullAssembly() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        //this->balance_ = eq_data_->balance_;
        this->element_cache_map_ = element_cache_map;

        shared_ptr<FE_P<dim>> fe_p = std::make_shared< FE_P<dim> >(1);
        shared_ptr<FE_P<dim-1>> fe_p_low = std::make_shared< FE_P<dim-1> >(1);
        fe_ = std::make_shared<FESystem<dim>>(fe_p, FEVector, 3);
        fe_low_ = std::make_shared<FESystem<dim-1>>(fe_p_low, FEVector, 3);
        fe_values_.initialize(*this->quad_, *fe_,
                update_values | update_gradients | update_JxW_values | update_quadrature_points);
        fe_values_side_.initialize(*this->quad_low_, *fe_,
                update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points);
        fe_values_sub_.initialize(*this->quad_low_, *fe_low_,
                update_values | update_gradients | update_JxW_values | update_quadrature_points);

        n_dofs_ = fe_->n_dofs();
        n_dofs_sub_ = fe_low_->n_dofs();
        n_dofs_ngh_ = { n_dofs_sub_, n_dofs_ };
        dof_indices_.resize(n_dofs_);
        side_dof_indices_.resize(2);
        side_dof_indices_[0].resize(n_dofs_sub_);  // index 0 = element with lower dimension,
        side_dof_indices_[1].resize(n_dofs_);      // index 1 = side of element with higher dimension
        local_matrix_.resize(n_dofs_*n_dofs_);
        local_matrix_ngh_.resize(2);
        for (uint m=0; m<2; ++m) {
            local_matrix_ngh_[m].resize(2);
            for (uint n=0; n<2; ++n)
                local_matrix_ngh_[m][n].resize(n_dofs_*n_dofs_);
        }
        vec_view_ = &fe_values_.vector_view(0);
        vec_view_side_ = &fe_values_side_.vector_view(0);
        if (dim>1) vec_view_sub_ = &fe_values_sub_.vector_view(0);
    }


    /// Assemble integral over element
    inline virtual void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        if (cell.dim() != dim) return;

        ElementAccessor<3> elm_acc = cell.elm();

        fe_values_.reinit(elm_acc);
        cell.get_dof_indices(dof_indices_);

        // assemble the local stiffness matrix
        for (unsigned int i=0; i<n_dofs_; i++)
            for (unsigned int j=0; j<n_dofs_; j++)
                local_matrix_[i*n_dofs_+j] = 0;

        unsigned int k=0;
        for (auto p : this->bulk_points(element_patch_idx) )
        {
            for (unsigned int i=0; i<n_dofs_; i++)
            {
                for (unsigned int j=0; j<n_dofs_; j++)
                    local_matrix_[i*n_dofs_+j] += eq_fields_->cross_section(p)*(
                                                2*eq_fields_->lame_mu(p)*arma::dot(vec_view_->sym_grad(j,k), vec_view_->sym_grad(i,k))
                                                + eq_fields_->lame_lambda(p)*vec_view_->divergence(j,k)*vec_view_->divergence(i,k)
                                               )*fe_values_.JxW(k);
            }
            k++;
        }
        this->cell_integral_set_values();
    }

    /// Assembles boundary integral.
    inline virtual void boundary_side_integral(DHCellSide cell_side)
    {
    	ASSERT_EQ(cell_side.dim(), dim).error("Dimension of element mismatch!");
        if (!cell_side.cell().is_own()) return;

        Side side = cell_side.side();
        const DHCellAccessor &dh_cell = cell_side.cell();
        dh_cell.get_dof_indices(dof_indices_);
        fe_values_side_.reinit(side);

        for (unsigned int i=0; i<n_dofs_; i++)
            for (unsigned int j=0; j<n_dofs_; j++)
                local_matrix_[i*n_dofs_+j] = 0;

        auto p_side = *( this->boundary_points(cell_side).begin() );
        auto p_bdr = p_side.point_bdr( side.cond().element_accessor() );
        unsigned int bc_type = eq_fields_->bc_type(p_bdr);
        double side_measure = cell_side.measure();
        if (bc_type == EqFields::bc_type_displacement)
        {
            unsigned int k=0;
            for (auto p : this->boundary_points(cell_side) ) {
                for (unsigned int i=0; i<n_dofs_; i++)
                    for (unsigned int j=0; j<n_dofs_; j++)
                        local_matrix_[i*n_dofs_+j] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
                                arma::dot(vec_view_side_->value(i,k),vec_view_side_->value(j,k)) * fe_values_side_.JxW(k);
                k++;
            }
        }
        else if (bc_type == EqFields::bc_type_displacement_normal)
        {
            unsigned int k=0;
            for (auto p : this->boundary_points(cell_side) ) {
                for (unsigned int i=0; i<n_dofs_; i++)
                    for (unsigned int j=0; j<n_dofs_; j++)
                        local_matrix_[i*n_dofs_+j] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
                                arma::dot(vec_view_side_->value(i,k), fe_values_side_.normal_vector(k)) *
                                arma::dot(vec_view_side_->value(j,k), fe_values_side_.normal_vector(k)) * fe_values_side_.JxW(k);
                k++;
            }
        }

        this->boundary_side_integral_set_values();
    }


    /// Assembles between elements of different dimensions.
    inline virtual void dimjoin_intergral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
    	if (dim == 1) return;
        ASSERT_EQ(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

		cell_lower_dim.get_dof_indices(side_dof_indices_[0]);
		ElementAccessor<3> cell_sub = cell_lower_dim.elm();
		fe_values_sub_.reinit(cell_sub);

		DHCellAccessor cell_higher_dim = eq_data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
		cell_higher_dim.get_dof_indices(side_dof_indices_[1]);
		fe_values_side_.reinit(neighb_side.side());

		// Element id's for testing if they belong to local partition.
		bool own_element_id[2];
		own_element_id[0] = cell_lower_dim.is_own();
		own_element_id[1] = cell_higher_dim.is_own();

        for (unsigned int n=0; n<2; ++n)
            for (unsigned int i=0; i<n_dofs_; i++)
                for (unsigned int m=0; m<2; ++m)
                    for (unsigned int j=0; j<n_dofs_; j++)
                        local_matrix_ngh_[n][m][i*(n_dofs_)+j] = 0;

        // set transmission conditions
        unsigned int k=0;
        for (auto p_high : this->coupling_points(neighb_side) )
        {
            auto p_low = p_high.lower_dim(cell_lower_dim);
            arma::vec3 nv = fe_values_side_.normal_vector(k);

            for (int n=0; n<2; n++)
            {
                if (!own_element_id[n]) continue;

                for (unsigned int i=0; i<n_dofs_ngh_[n]; i++)
                {
                    arma::vec3 vi = (n==0) ? arma::zeros(3) : vec_view_side_->value(i,k);
                    arma::vec3 vf = (n==1) ? arma::zeros(3) : vec_view_sub_->value(i,k);
                    arma::mat33 gvft = (n==0) ? vec_view_sub_->grad(i,k) : arma::zeros(3,3);

                    for (int m=0; m<2; m++)
                        for (unsigned int j=0; j<n_dofs_ngh_[m]; j++) {
                            arma::vec3 ui = (m==0) ? arma::zeros(3) : vec_view_side_->value(j,k);
                            arma::vec3 uf = (m==1) ? arma::zeros(3) : vec_view_sub_->value(j,k);
                            arma::mat33 guft = (m==0) ? mat_t(vec_view_sub_->grad(j,k),nv) : arma::zeros(3,3);
                            double divuft = (m==0) ? arma::trace(guft) : 0;

                            local_matrix_ngh_[n][m][i*n_dofs_ngh_[m] + j] +=
                                    eq_fields_->fracture_sigma(p_low)*(
                                     arma::dot(vf-vi,
                                      2/eq_fields_->cross_section(p_low)*(eq_fields_->lame_mu(p_low)*(uf-ui)+(eq_fields_->lame_mu(p_low)+eq_fields_->lame_lambda(p_low))*(arma::dot(uf-ui,nv)*nv))
                                      + eq_fields_->lame_mu(p_low)*arma::trans(guft)*nv
                                      + eq_fields_->lame_lambda(p_low)*divuft*nv
                                     )
                                     - arma::dot(gvft, eq_fields_->lame_mu(p_low)*arma::kron(nv,ui.t()) + eq_fields_->lame_lambda(p_low)*arma::dot(ui,nv)*arma::eye(3,3))
                                    )*fe_values_sub_.JxW(k);
                        }

                }
            }
        	k++;
        }

        this->dimjoin_intergral_set_values();
    }



protected:
    inline virtual void cell_integral_set_values() {
        this->eq_data_->ls->mat_set_values(n_dofs_, dof_indices_.data(), n_dofs_, dof_indices_.data(), &(local_matrix_[0]));
    }

    inline virtual void boundary_side_integral_set_values() {
        this->eq_data_->ls->mat_set_values(n_dofs_, dof_indices_.data(), n_dofs_, dof_indices_.data(), &(local_matrix_[0]));
    }

    inline virtual void dimjoin_intergral_set_values() {
        for (unsigned int n=0; n<2; ++n)
            for (unsigned int m=0; m<2; ++m)
                this->eq_data_->ls->mat_set_values(n_dofs_ngh_[n], side_dof_indices_[n].data(), n_dofs_ngh_[m], side_dof_indices_[m].data(), &(local_matrix_ngh_[n][m][0]));
    }


    inline arma::mat33 mat_t(const arma::mat33 &m, const arma::vec3 &n)
    {
      arma::mat33 mt = m - m*arma::kron(n,n.t());
      return mt;
    }



    shared_ptr<FiniteElement<dim>> fe_;         ///< Finite element for the solution of the advection-diffusion equation.
    shared_ptr<FiniteElement<dim-1>> fe_low_;   ///< Finite element for the solution of the advection-diffusion equation (dim-1).

    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int n_dofs_;                                     ///< Number of dofs
    unsigned int n_dofs_sub_;                                 ///< Number of dofs (on lower dim element)
    std::vector<unsigned int> n_dofs_ngh_;                    ///< Number of dofs on lower and higher dimension element (vector of 2 items)
    FEValues<3> fe_values_;                                   ///< FEValues of cell object (FESystem of P disc finite element type)
    FEValues<3> fe_values_side_;                              ///< FEValues of side object
    FEValues<3> fe_values_sub_;                               ///< FEValues of lower dimension cell object

    vector<LongIdx> dof_indices_;                             ///< Vector of global DOF indices
    vector<vector<LongIdx> > side_dof_indices_;               ///< 2 items vector of DOF indices in neighbour calculation.
    vector<PetscScalar> local_matrix_;                        ///< Auxiliary vector for assemble methods
    vector<vector<vector<PetscScalar>>> local_matrix_ngh_;    ///< Auxiliary vectors for assemble ngh integral
    const FEValuesViews::Vector<3> * vec_view_;               ///< Vector view in cell integral calculation.
    const FEValuesViews::Vector<3> * vec_view_side_;          ///< Vector view in boundary / neighbour calculation.
    const FEValuesViews::Vector<3> * vec_view_sub_;           ///< Vector view of low dim element in neighbour calculation.

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


template <unsigned int dim>
class Stiffness_ComputeLocal : public Stiffness_FullAssembly<dim>
{
public:
    typedef equation_data::EqFields EqFields;
    typedef equation_data::EqData EqData;

    static constexpr const char * name() { return "StiffnessAssembly"; }

    /// Constructor.
    Stiffness_ComputeLocal(EqFields *eq_fields, EqData *eq_data)
    : Stiffness_FullAssembly<dim>(eq_fields, eq_data) {}

    /// Destructor.
    ~Stiffness_ComputeLocal() {}

protected:
    void cell_integral_set_values() override {}

    void boundary_side_integral_set_values() override {}

    void dimjoin_intergral_set_values() override {}

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


template <unsigned int dim>
class Stiffness_EvalFields : public Stiffness_FullAssembly<dim>
{
public:
    typedef equation_data::EqFields EqFields;
    typedef equation_data::EqData EqData;

    static constexpr const char * name() { return "StiffnessAssembly"; }

    /// Constructor.
    Stiffness_EvalFields(EqFields *eq_fields, EqData *eq_data)
    : Stiffness_FullAssembly<dim>(eq_fields, eq_data) {}

    /// Destructor.
    ~Stiffness_EvalFields() {}

    /// Assembles the cell (volume) integral into the stiffness matrix.
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx) override {}

    /// Assembles the fluxes on the boundary.
    inline void boundary_side_integral(DHCellSide cell_side) override {}

    /// Assembles the fluxes between elements of different dimensions.
    inline void dimjoin_intergral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) override {}

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


template <unsigned int dim>
class Rhs_FullAssembly : public AssemblyBase<dim>
{
public:
    typedef typename equation_data::EqFields EqFields;
    typedef typename equation_data::EqData EqData;

    static constexpr const char * name() { return "RhsAssembly"; }

    /// Constructor.
    Rhs_FullAssembly(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(1), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = (ActiveIntegrals::bulk | ActiveIntegrals::coupling | ActiveIntegrals::boundary);
        this->used_fields_ += eq_fields_->cross_section;
        this->used_fields_ += eq_fields_->load;
        this->used_fields_ += eq_fields_->potential_load;
        this->used_fields_ += eq_fields_->ref_potential_load;
        this->used_fields_ += eq_fields_->fracture_sigma;
        this->used_fields_ += eq_fields_->dirichlet_penalty;
        this->used_fields_ += eq_fields_->bc_type;
        this->used_fields_ += eq_fields_->bc_displacement;
        this->used_fields_ += eq_fields_->bc_traction;
        this->used_fields_ += eq_fields_->bc_stress;
        this->used_fields_ += eq_fields_->initial_stress;
    }

    /// Destructor.
    ~Rhs_FullAssembly() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        //this->balance_ = eq_data_->balance_;
        this->element_cache_map_ = element_cache_map;

        shared_ptr<FE_P<dim>> fe_p = std::make_shared< FE_P<dim> >(1);
        shared_ptr<FE_P<dim-1>> fe_p_low = std::make_shared< FE_P<dim-1> >(1);
        fe_ = std::make_shared<FESystem<dim>>(fe_p, FEVector, 3);
        fe_low_ = std::make_shared<FESystem<dim-1>>(fe_p_low, FEVector, 3);
        fe_values_.initialize(*this->quad_, *fe_,
                update_values | update_gradients | update_JxW_values | update_quadrature_points);
        fe_values_bdr_side_.initialize(*this->quad_low_, *fe_,
                update_values | update_normal_vectors | update_side_JxW_values | update_quadrature_points);
        fe_values_side_.initialize(*this->quad_low_, *fe_,
                update_values | update_normal_vectors);
        fe_values_sub_.initialize(*this->quad_low_, *fe_low_,
                update_values | update_JxW_values | update_quadrature_points);
        n_dofs_ = fe_->n_dofs();
        n_dofs_sub_ = fe_low_->n_dofs();
        n_dofs_ngh_ = { n_dofs_sub_, n_dofs_ };
        dof_indices_.resize(n_dofs_);
        side_dof_indices_.resize(2);
        side_dof_indices_[0].resize(n_dofs_sub_);  // index 0 = element with lower dimension,
        side_dof_indices_[1].resize(n_dofs_);      // index 1 = side of element with higher dimension
        local_rhs_.resize(n_dofs_);
        local_rhs_ngh_.resize(2);
        for (uint n=0; n<2; ++n) local_rhs_ngh_[n].resize(n_dofs_);
        vec_view_ = &fe_values_.vector_view(0);
        vec_view_bdr_ = &fe_values_bdr_side_.vector_view(0);
        vec_view_side_ = &fe_values_side_.vector_view(0);
        if (dim>1) vec_view_sub_ = &fe_values_sub_.vector_view(0);
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        if (cell.dim() != dim) return;
        if (!cell.is_own()) return;

        ElementAccessor<3> elm_acc = cell.elm();

        fe_values_.reinit(elm_acc);
        cell.get_dof_indices(dof_indices_);

        // assemble the local stiffness matrix
        fill_n(&(local_rhs_[0]), n_dofs_, 0);

        // compute sources
        unsigned int k=0;
        for (auto p : this->bulk_points(element_patch_idx) )
        {
            for (unsigned int i=0; i<n_dofs_; i++)
                local_rhs_[i] += (
                                 arma::dot(eq_fields_->load(p), vec_view_->value(i,k))
                                 -eq_fields_->potential_load(p)*vec_view_->divergence(i,k)
                                 -arma::dot(eq_fields_->initial_stress(p), vec_view_->grad(i,k))
                                )*eq_fields_->cross_section(p)*fe_values_.JxW(k);
            ++k;
        }
        this->cell_integral_set_values();
    }

    /// Assembles boundary integral.
    inline void boundary_side_integral(DHCellSide cell_side)
    {
    	ASSERT_EQ(cell_side.dim(), dim).error("Dimension of element mismatch!");
        if (!cell_side.cell().is_own()) return;

        Side side = cell_side.side();
        const DHCellAccessor &dh_cell = cell_side.cell();
        dh_cell.get_dof_indices(dof_indices_);
        fe_values_bdr_side_.reinit(side);

        auto p_side = *( this->boundary_points(cell_side).begin() );
        auto p_bdr = p_side.point_bdr( side.cond().element_accessor() );
        unsigned int bc_type = eq_fields_->bc_type(p_bdr);

        fill_n(&(local_rhs_[0]), n_dofs_, 0);

        unsigned int k = 0;

        // addtion from initial stress
        for (auto p : this->boundary_points(cell_side) )
        {
            for (unsigned int i=0; i<n_dofs_; i++)
                local_rhs_[i] += eq_fields_->cross_section(p) *
                        arma::dot(( eq_fields_->initial_stress(p) * fe_values_bdr_side_.normal_vector(k)),
                                    vec_view_bdr_->value(i,k)) *
                        fe_values_bdr_side_.JxW(k);
            ++k;
        }

        k = 0;
        if (bc_type == EqFields::bc_type_displacement)
        {
            double side_measure = cell_side.measure();
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( side.cond().element_accessor() );
                for (unsigned int i=0; i<n_dofs_; i++)
                    local_rhs_[i] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
					        arma::dot(eq_fields_->bc_displacement(p_bdr), vec_view_bdr_->value(i,k)) *
					        fe_values_bdr_side_.JxW(k);
                ++k;
            }
        }
        else if (bc_type == EqFields::bc_type_displacement_normal)
        {
            double side_measure = cell_side.measure();
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( side.cond().element_accessor() );
                for (unsigned int i=0; i<n_dofs_; i++)
                    local_rhs_[i] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
                            arma::dot(eq_fields_->bc_displacement(p_bdr), fe_values_bdr_side_.normal_vector(k)) *
                            arma::dot(vec_view_bdr_->value(i,k), fe_values_bdr_side_.normal_vector(k)) *
                            fe_values_bdr_side_.JxW(k);
                ++k;
            }
        }
        else if (bc_type == EqFields::bc_type_traction)
        {
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( side.cond().element_accessor() );
                for (unsigned int i=0; i<n_dofs_; i++)
                    local_rhs_[i] += eq_fields_->cross_section(p) *
                            arma::dot(vec_view_bdr_->value(i,k), eq_fields_->bc_traction(p_bdr) + eq_fields_->ref_potential_load(p) * fe_values_bdr_side_.normal_vector(k)) *
                            fe_values_bdr_side_.JxW(k);
                ++k;
            }
        }
        else if (bc_type == EqFields::bc_type_stress)
        {
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( side.cond().element_accessor() );
                for (unsigned int i=0; i<n_dofs_; i++)
                    // stress is multiplied by inward normal to obtain traction
                    local_rhs_[i] += eq_fields_->cross_section(p) *
                            arma::dot(vec_view_bdr_->value(i,k), -eq_fields_->bc_stress(p_bdr)*fe_values_bdr_side_.normal_vector(k)
                            + eq_fields_->ref_potential_load(p) * fe_values_bdr_side_.normal_vector(k))
                            * fe_values_bdr_side_.JxW(k);
                ++k;
            }
        }
        this->boundary_side_integral_set_values();
    }


    /// Assembles between elements of different dimensions.
    inline void dimjoin_intergral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
    	if (dim == 1) return;
        ASSERT_EQ(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

		cell_lower_dim.get_dof_indices(side_dof_indices_[0]);
		ElementAccessor<3> cell_sub = cell_lower_dim.elm();
		fe_values_sub_.reinit(cell_sub);

		DHCellAccessor cell_higher_dim = eq_data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
		cell_higher_dim.get_dof_indices(side_dof_indices_[1]);
		fe_values_side_.reinit(neighb_side.side());

		// Element id's for testing if they belong to local partition.
		bool own_element_id[2];
		own_element_id[0] = cell_lower_dim.is_own();
		own_element_id[1] = cell_higher_dim.is_own();

        for (unsigned int n=0; n<2; ++n)
            for (unsigned int i=0; i<n_dofs_; i++)
                local_rhs_ngh_[n][i] = 0;

        // set transmission conditions
        unsigned int k=0;
        for (auto p_high : this->coupling_points(neighb_side) )
        {
            auto p_low = p_high.lower_dim(cell_lower_dim);
            arma::vec3 nv = fe_values_side_.normal_vector(k);

            for (int n=0; n<2; n++)
            {
                if (!own_element_id[n]) continue;

                for (unsigned int i=0; i<n_dofs_ngh_[n]; i++)
                {
                    arma::vec3 vi = (n==0) ? arma::zeros(3) : vec_view_side_->value(i,k);
                    arma::vec3 vf = (n==1) ? arma::zeros(3) : vec_view_sub_->value(i,k);

                    local_rhs_ngh_[n][i] -= eq_fields_->fracture_sigma(p_low) * eq_fields_->cross_section(p_high) *
                            arma::dot(vf-vi, eq_fields_->potential_load(p_high) * nv) * fe_values_sub_.JxW(k);
                }
            }
            ++k;
        }

        this->dimjoin_intergral_set_values();
    }



protected:
    inline virtual void cell_integral_set_values() {
        this->eq_data_->ls->rhs_set_values(n_dofs_, dof_indices_.data(), &(local_rhs_[0]));
    }

    inline virtual void boundary_side_integral_set_values() {
        this->eq_data_->ls->rhs_set_values(n_dofs_, dof_indices_.data(), &(local_rhs_[0]));
    }

    inline virtual void dimjoin_intergral_set_values() {
        for (unsigned int n=0; n<2; ++n)
            this->eq_data_->ls->rhs_set_values(n_dofs_ngh_[n], side_dof_indices_[n].data(), &(local_rhs_ngh_[n][0]));
    }


    shared_ptr<FiniteElement<dim>> fe_;         ///< Finite element for the solution of the advection-diffusion equation.
    shared_ptr<FiniteElement<dim-1>> fe_low_;   ///< Finite element for the solution of the advection-diffusion equation (dim-1).

    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int n_dofs_;                                     ///< Number of dofs
    unsigned int n_dofs_sub_;                                 ///< Number of dofs (on lower dim element)
    std::vector<unsigned int> n_dofs_ngh_;                    ///< Number of dofs on lower and higher dimension element (vector of 2 items)
    FEValues<3> fe_values_;                                   ///< FEValues of cell object (FESystem of P disc finite element type)
    FEValues<3> fe_values_bdr_side_;                          ///< FEValues of side (boundary integral) object
    FEValues<3> fe_values_side_;                              ///< FEValues of side (neighbour integral) object
    FEValues<3> fe_values_sub_;                               ///< FEValues of lower dimension cell object

    vector<LongIdx> dof_indices_;                             ///< Vector of global DOF indices
    vector<vector<LongIdx> > side_dof_indices_;               ///< 2 items vector of DOF indices in neighbour calculation.
    vector<PetscScalar> local_rhs_;                           ///< Auxiliary vector for assemble methods
    vector<vector<PetscScalar>> local_rhs_ngh_;               ///< Auxiliary vectors for assemble ngh integral
    const FEValuesViews::Vector<3> * vec_view_;               ///< Vector view in cell integral calculation.
    const FEValuesViews::Vector<3> * vec_view_bdr_;           ///< Vector view in boundary calculation.
    const FEValuesViews::Vector<3> * vec_view_side_;          ///< Vector view in neighbour calculation.
    const FEValuesViews::Vector<3> * vec_view_sub_;           ///< Vector view of low dim element in neighbour calculation.


    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


template <unsigned int dim>
class Rhs_ComputeLocal : public Rhs_FullAssembly<dim>
{
public:
    typedef equation_data::EqFields EqFields;
    typedef equation_data::EqData EqData;

    static constexpr const char * name() { return "RhsAssembly"; }

    /// Constructor.
    Rhs_ComputeLocal(EqFields *eq_fields, EqData *eq_data)
    : Rhs_FullAssembly<dim>(eq_fields, eq_data) {}

    /// Destructor.
    ~Rhs_ComputeLocal() {}

protected:
    void cell_integral_set_values() override {}

    void boundary_side_integral_set_values() override {}

    void dimjoin_intergral_set_values() override {}

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


template <unsigned int dim>
class Rhs_EvalFields : public Rhs_FullAssembly<dim>
{
public:
    typedef equation_data::EqFields EqFields;
    typedef equation_data::EqData EqData;

    static constexpr const char * name() { return "RhsAssembly"; }

    /// Constructor.
    Rhs_EvalFields(EqFields *eq_fields, EqData *eq_data)
    : Rhs_FullAssembly<dim>(eq_fields, eq_data) {}

    /// Destructor.
    ~Rhs_EvalFields() {}

    /// Assembles the cell (volume) integral into the stiffness matrix.
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx) override {}

    /// Assembles the fluxes on the boundary.
    inline void boundary_side_integral(DHCellSide cell_side) override {}

    /// Assembles the fluxes between elements of different dimensions.
    inline void dimjoin_intergral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) override {}

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


#endif /* ELASTICITY_MOCKUP_ASSEMBLYA_HH_ */

