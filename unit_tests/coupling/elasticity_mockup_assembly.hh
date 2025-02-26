#ifndef ELASTICITY_MOCKUP_ASSEMBLYA_HH_
#define ELASTICITY_MOCKUP_ASSEMBLYA_HH_

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "elasticity_mockup.hh"
#include "fem/fe_p.hh"
#include "fem/patch_fe_values.hh"
#include "quadrature/quadrature_lib.hh"
#include "coupling/balance.hh"
#include "fem/element_cache_map.hh"


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim>
class Stiffness_FullAssembly : public AssemblyBasePatch<dim>
{
public:
    typedef typename equation_data::EqFields EqFields;
    typedef typename equation_data::EqData EqData;

    static constexpr const char * name() { return "StiffnessAssembly"; }

    /// Constructor.
    Stiffness_FullAssembly(EqFields *eq_fields, EqData *eq_data, PatchFEValues<3> *fe_values)
    : AssemblyBasePatch<dim>(fe_values), eq_fields_(eq_fields), eq_data_(eq_data), // quad_order = 1
      JxW_( this->bulk_values().JxW() ),
      JxW_side_( this->side_values().JxW() ),
      normal_( this->side_values().normal_vector() ),
      deform_side_( this->side_values().vector_shape() ),
      grad_deform_( this->bulk_values().grad_vector_shape() ),
      sym_grad_deform_( this->bulk_values().vector_sym_grad() ),
      div_deform_( this->bulk_values().vector_divergence() ),
      deform_join_( this->join_values().vector_join_shape() ),
      deform_join_grad_( this->join_values().gradient_vector_join_shape() ) {
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

        shared_ptr<FE_P<dim-1>> fe_p_low = std::make_shared< FE_P<dim-1> >(1);
        shared_ptr<FiniteElement<dim-1>> fe_low = std::make_shared<FESystem<dim-1>>(fe_p_low, FEVector, 3);
        this->fe_values_->template initialize<dim>(*this->quad_);
        this->fe_values_->template initialize<dim>(*this->quad_low_);

        n_dofs_ = this->n_dofs();
        n_dofs_sub_ = fe_low->n_dofs();
        n_dofs_ngh_ = { n_dofs_sub_, n_dofs_ };
        dof_indices_.resize(n_dofs_);
        side_dof_indices_.resize(2*n_dofs_);
        local_matrix_.resize(4*n_dofs_*n_dofs_);
    }


    /// Assemble integral over element
    inline virtual void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        if (cell.dim() != dim) return;

        // Fracture stiffness is assembled in dimjoin_integral.
        if (cell.elm()->n_neighs_vb() > 0) return;

        cell.get_dof_indices(dof_indices_);

        // assemble the local stiffness matrix
        for (unsigned int i=0; i<n_dofs_; i++)
            for (unsigned int j=0; j<n_dofs_; j++)
                local_matrix_[i*n_dofs_+j] = 0;

        for (auto p : this->bulk_points(element_patch_idx) )
        {
            for (unsigned int i=0; i<n_dofs_; i++)
            {
                for (unsigned int j=0; j<n_dofs_; j++)
                    local_matrix_[i*n_dofs_+j] += eq_fields_->cross_section(p)*(
                                                arma::dot(eq_fields_->stress_tensor(p,sym_grad_deform_.shape(j)(p)), sym_grad_deform_.shape(i)(p))
                                               )*JxW_(p);
            }
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

        for (unsigned int i=0; i<n_dofs_; i++)
            for (unsigned int j=0; j<n_dofs_; j++)
                local_matrix_[i*n_dofs_+j] = 0;

        auto p_side = *( this->boundary_points(cell_side).begin() );
        auto p_bdr = p_side.point_bdr( side.cond().element_accessor() );
        unsigned int bc_type = eq_fields_->bc_type(p_bdr);
        double side_measure = cell_side.measure();
        if (bc_type == EqFields::bc_type_displacement)
        {
            for (auto p : this->boundary_points(cell_side) ) {
                for (unsigned int i=0; i<n_dofs_; i++)
                    for (unsigned int j=0; j<n_dofs_; j++)
                        local_matrix_[i*n_dofs_+j] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
                                arma::dot(deform_side_.shape(i)(p),deform_side_.shape(j)(p)) * JxW_side_(p);
            }
        }
        else if (bc_type == EqFields::bc_type_displacement_normal)
        {
            for (auto p : this->boundary_points(cell_side) ) {
                for (unsigned int i=0; i<n_dofs_; i++)
                    for (unsigned int j=0; j<n_dofs_; j++)
                        local_matrix_[i*n_dofs_+j] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
                                arma::dot(deform_side_.shape(i)(p), normal_(p)) *
                                arma::dot(deform_side_.shape(j)(p), normal_(p)) * JxW_side_(p);
            }
        }

        this->boundary_side_integral_set_values();
    }


    /// Assembles between elements of different dimensions.
    inline virtual void dimjoin_intergral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
        if (dim == 1) return;
        ASSERT_EQ(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

        unsigned int n_indices = cell_lower_dim.get_dof_indices(dof_indices_);
        for(unsigned int i=0; i<n_indices; ++i) {
            side_dof_indices_[i] = dof_indices_[i];
        }

        DHCellAccessor cell_higher_dim = eq_data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
        n_indices = cell_higher_dim.get_dof_indices(dof_indices_);
        for(unsigned int i=0; i<n_indices; ++i) {
            side_dof_indices_[i+n_dofs_ngh_[0]] = dof_indices_[i];
        }

        // Element id's for testing if they belong to local partition.
        bool own_element_id[2];
        own_element_id[0] = cell_lower_dim.is_own();
        own_element_id[1] = cell_higher_dim.is_own();

        unsigned int n_neighs = cell_lower_dim.elm()->n_neighs_vb();

        for (unsigned int i=0; i<n_dofs_ngh_[0]+n_dofs_ngh_[1]; i++)
            for (unsigned int j=0; j<n_dofs_ngh_[0]+n_dofs_ngh_[1]; j++)
                local_matrix_[i*(n_dofs_ngh_[0]+n_dofs_ngh_[1])+j] = 0;

        // set transmission conditions
        for (auto p_high : this->coupling_points(neighb_side) )
        {
            auto p_low = p_high.lower_dim(cell_lower_dim);
            arma::vec3 nv = normal_(p_high);

            for (uint i=0; i<deform_join_.n_dofs_both(); ++i) {
                uint is_high_i = deform_join_.is_high_dim(i);
                if (!own_element_id[is_high_i]) continue;
                arma::vec3 diff_deform_i = deform_join_.shape(i)(p_low) - deform_join_.shape(i)(p_high);
                arma::mat33 grad_deform_i = deform_join_grad_.shape(i)(p_low);  // low dim element
                arma::mat33 semi_grad_i = grad_deform_i + n_neighs/eq_fields_->cross_section(p_low)*arma::kron(diff_deform_i,nv.t());
                arma::mat33 semi_sym_grad_i = 0.5*(semi_grad_i + semi_grad_i.t());

                for (uint j=0; j<deform_join_.n_dofs_both(); ++j) {
                    arma::vec3 deform_j_high = deform_join_.shape(j)(p_high);
                    arma::vec3 diff_deform_j = deform_join_.shape(j)(p_low) - deform_j_high;
                    arma::mat33 grad_deform_j = deform_join_grad_.shape(j)(p_low);  // low dim element
                    arma::mat33 semi_grad_j = grad_deform_j + n_neighs/eq_fields_->cross_section(p_low)*arma::kron(diff_deform_j,nv.t());
                    arma::mat33 semi_sym_grad_j = 0.5*(semi_grad_j + semi_grad_j.t());

                    local_matrix_[i * (n_dofs_ngh_[0]+n_dofs_ngh_[1]) + j] +=
                            (
                                     eq_fields_->fracture_sigma(p_low)*eq_fields_->cross_section(p_low) / n_neighs
                                     * arma::dot(semi_sym_grad_i, eq_fields_->stress_tensor(p_low,semi_sym_grad_j))

                                     // This term corrects the tangential part of the above product so that there is no
                                     // dependence on fracture_sigma.
                                     // TODO: Fracture_sigma should be possibly removed and replaced by anisotropic elasticity.
                                     + (1-eq_fields_->fracture_sigma(p_low))*eq_fields_->cross_section(p_low) / n_neighs
                                     * arma::dot(0.5*(grad_deform_i+grad_deform_i.t()), eq_fields_->stress_tensor(p_low,0.5*(grad_deform_j+grad_deform_j.t())))
                            )*JxW_side_(p_high);
                }
            }
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
        this->eq_data_->ls->mat_set_values(n_dofs_ngh_[0]+n_dofs_ngh_[1], &(side_dof_indices_[0]), n_dofs_ngh_[0]+n_dofs_ngh_[1], &(side_dof_indices_[0]), &(local_matrix_[0]));
        this->eq_data_->ls->mat_set_values(n_dofs_ngh_[0]+n_dofs_ngh_[1], &(side_dof_indices_[0]), n_dofs_ngh_[0]+n_dofs_ngh_[1], &(side_dof_indices_[0]), &(local_matrix_[0]));
    }


    inline arma::mat33 mat_t(const arma::mat33 &m, const arma::vec3 &n)
    {
      arma::mat33 mt = m - m*arma::kron(n,n.t());
      return mt;
    }



    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int n_dofs_;                                     ///< Number of dofs
    unsigned int n_dofs_sub_;                                 ///< Number of dofs (on lower dim element)
    std::vector<unsigned int> n_dofs_ngh_;                    ///< Number of dofs on lower and higher dimension element (vector of 2 items)

    vector<LongIdx> dof_indices_;                             ///< Vector of global DOF indices
    vector<LongIdx> side_dof_indices_;                        ///< 2 items vector of DOF indices in neighbour calculation.
    vector<PetscScalar> local_matrix_;                        ///< Auxiliary vector for assemble methods

    /// Following data members represent Element quantities and FE quantities
    FeQ<Scalar> JxW_;
    FeQ<Scalar> JxW_side_;
    ElQ<Vector> normal_;
    FeQArray<Vector> deform_side_;
    FeQArray<Tensor> grad_deform_;
    FeQArray<Tensor> sym_grad_deform_;
    FeQArray<Scalar> div_deform_;
    FeQJoin<Vector> deform_join_;
    FeQJoin<Tensor> deform_join_grad_;

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
    Stiffness_ComputeLocal(EqFields *eq_fields, EqData *eq_data, PatchFEValues<3> *fe_values)
    : Stiffness_FullAssembly<dim>(eq_fields, eq_data, fe_values) {}

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
    Stiffness_EvalFields(EqFields *eq_fields, EqData *eq_data, PatchFEValues<3> *fe_values)
    : Stiffness_FullAssembly<dim>(eq_fields, eq_data, fe_values) {}

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
class Rhs_FullAssembly : public AssemblyBasePatch<dim>
{
public:
    typedef typename equation_data::EqFields EqFields;
    typedef typename equation_data::EqData EqData;

    static constexpr const char * name() { return "RhsAssembly"; }

    /// Constructor.
    Rhs_FullAssembly(EqFields *eq_fields, EqData *eq_data, PatchFEValues<3> *fe_values)
    : AssemblyBasePatch<dim>(fe_values), eq_fields_(eq_fields), eq_data_(eq_data),
      JxW_( this->bulk_values().JxW() ),
      JxW_side_( this->side_values().JxW() ),
      normal_( this->side_values().normal_vector() ),
      deform_( this->bulk_values().vector_shape() ),
      deform_side_( this->side_values().vector_shape() ),
	  grad_deform_( this->bulk_values().grad_vector_shape() ),
      div_deform_( this->bulk_values().vector_divergence() ),
      deform_join_( this->join_values().vector_join_shape() ) {
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

        shared_ptr<FE_P<dim-1>> fe_p_low = std::make_shared< FE_P<dim-1> >(1);
        shared_ptr<FiniteElement<dim-1>> fe_low = std::make_shared<FESystem<dim-1>>(fe_p_low, FEVector, 3);
        this->fe_values_->template initialize<dim>(*this->quad_);
        this->fe_values_->template initialize<dim>(*this->quad_low_);

        n_dofs_ = this->n_dofs();
        n_dofs_sub_ = fe_low->n_dofs();
        n_dofs_ngh_ = { n_dofs_sub_, n_dofs_ };
        dof_indices_.resize(n_dofs_);
        side_dof_indices_.resize(n_dofs_sub_ + n_dofs_);
        local_rhs_.resize(2*n_dofs_);
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        if (cell.dim() != dim) return;
        if (!cell.is_own()) return;

        cell.get_dof_indices(dof_indices_);

        // assemble the local stiffness matrix
        fill_n(&(local_rhs_[0]), n_dofs_, 0);
        //local_source_balance_vector.assign(n_dofs_, 0);
        //local_source_balance_rhs.assign(n_dofs_, 0);

        // compute sources
        for (auto p : this->bulk_points(element_patch_idx) )
        {
            for (unsigned int i=0; i<n_dofs_; i++)
                local_rhs_[i] += (
                                 arma::dot(eq_fields_->load(p), deform_.shape(i)(p))
                                 -eq_fields_->potential_load(p)*div_deform_.shape(i)(p)
                                 -arma::dot(eq_fields_->initial_stress(p), grad_deform_.shape(i)(p))
                                )*eq_fields_->cross_section(p)*JxW_(p);
        }
        this->cell_integral_set_values();
    }

    /// Assembles boundary integral.
    inline void boundary_side_integral(DHCellSide cell_side)
    {
    	ASSERT_EQ(cell_side.dim(), dim).error("Dimension of element mismatch!");
        if (!cell_side.cell().is_own()) return;

        const DHCellAccessor &dh_cell = cell_side.cell();
        dh_cell.get_dof_indices(dof_indices_);

        auto p_side = *( this->boundary_points(cell_side).begin() );
        auto p_bdr = p_side.point_bdr( cell_side.cond().element_accessor() );
        unsigned int bc_type = eq_fields_->bc_type(p_bdr);

        fill_n(&(local_rhs_[0]), n_dofs_, 0);
        // local_flux_balance_vector.assign(n_dofs_, 0);
        // local_flux_balance_rhs = 0;

        // addtion from initial stress
        for (auto p : this->boundary_points(cell_side) )
        {
            for (unsigned int i=0; i<n_dofs_; i++)
                local_rhs_[i] += eq_fields_->cross_section(p) *
                        arma::dot(( eq_fields_->initial_stress(p) * normal_(p)),
                                    deform_side_.shape(i)(p)) *
                        JxW_side_(p);
        }

        if (bc_type == EqFields::bc_type_displacement)
        {
            double side_measure = cell_side.measure();
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( cell_side.cond().element_accessor() );
                for (unsigned int i=0; i<n_dofs_; i++)
                    local_rhs_[i] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
					        arma::dot(eq_fields_->bc_displacement(p_bdr), deform_side_.shape(i)(p)) *
					        JxW_side_(p);
            }
        }
        else if (bc_type == EqFields::bc_type_displacement_normal)
        {
            double side_measure = cell_side.measure();
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( cell_side.cond().element_accessor() );
                for (unsigned int i=0; i<n_dofs_; i++)
                    local_rhs_[i] += (eq_fields_->dirichlet_penalty(p) / side_measure) *
                            arma::dot(eq_fields_->bc_displacement(p_bdr), normal_(p)) *
                            arma::dot(deform_side_.shape(i)(p), normal_(p)) *
                            JxW_side_(p);
            }
        }
        else if (bc_type == EqFields::bc_type_traction)
        {
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( cell_side.cond().element_accessor() );
                for (unsigned int i=0; i<n_dofs_; i++)
                    local_rhs_[i] += eq_fields_->cross_section(p) *
                            arma::dot(deform_side_.shape(i)(p), eq_fields_->bc_traction(p_bdr) + eq_fields_->ref_potential_load(p) * normal_(p)) *
                            JxW_side_(p);
            }
        }
        else if (bc_type == EqFields::bc_type_stress)
        {
            for (auto p : this->boundary_points(cell_side) )
            {
                auto p_bdr = p.point_bdr( cell_side.cond().element_accessor() );
                for (unsigned int i=0; i<n_dofs_; i++)
                    // stress is multiplied by inward normal to obtain traction
                    local_rhs_[i] += eq_fields_->cross_section(p) *
                            arma::dot(deform_side_.shape(i)(p), -eq_fields_->bc_stress(p_bdr)*normal_(p)
                            + eq_fields_->ref_potential_load(p) * normal_(p))
                            * JxW_side_(p);
            }
        }
        this->boundary_side_integral_set_values();
    }


    /// Assembles between elements of different dimensions.
    inline void dimjoin_intergral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
    	if (dim == 1) return;
        ASSERT_EQ(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

        unsigned int n_indices = cell_lower_dim.get_dof_indices(dof_indices_);
        for(unsigned int i=0; i<n_indices; ++i) {
            side_dof_indices_[i] = dof_indices_[i];
        }

	    DHCellAccessor cell_higher_dim = eq_data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
        n_indices = cell_higher_dim.get_dof_indices(dof_indices_);
        for(unsigned int i=0; i<n_indices; ++i) {
            side_dof_indices_[i+n_dofs_ngh_[0]] = dof_indices_[i];
        }

		// Element id's for testing if they belong to local partition.
		bool own_element_id[2];
		own_element_id[0] = cell_lower_dim.is_own();
		own_element_id[1] = cell_higher_dim.is_own();

        for (unsigned int i=0; i<2*n_dofs_; i++)
            local_rhs_[i] = 0;

        // set transmission conditions
        for (auto p_high : this->coupling_points(neighb_side) )
        {
            auto p_low = p_high.lower_dim(cell_lower_dim);
            arma::vec3 nv = normal_(p_high);

            for (uint i=0; i<deform_join_.n_dofs_both(); ++i) {
                uint is_high_i = deform_join_.is_high_dim(i);
                if (!own_element_id[is_high_i]) continue;

                arma::vec3 vi = deform_join_.shape(i)(p_high);
                arma::vec3 vf = deform_join_.shape(i)(p_low);

                local_rhs_[i] -= eq_fields_->fracture_sigma(p_low) * eq_fields_->cross_section(p_high) *
                        arma::dot(vf-vi, eq_fields_->potential_load(p_high) * nv) * JxW_side_(p_high);
            }
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
        this->eq_data_->ls->rhs_set_values(n_dofs_ngh_[0]+n_dofs_ngh_[1], side_dof_indices_.data(), &(local_rhs_[0]));
    }


    /// Data objects shared with Elasticity
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int n_dofs_;                                     ///< Number of dofs
    unsigned int n_dofs_sub_;                                 ///< Number of dofs (on lower dim element)
    std::vector<unsigned int> n_dofs_ngh_;                    ///< Number of dofs on lower and higher dimension element (vector of 2 items)

    vector<LongIdx> dof_indices_;                                       ///< Vector of global DOF indices
    vector<LongIdx> side_dof_indices_;                                  ///< 2 items vector of DOF indices in neighbour calculation.
    vector<PetscScalar> local_rhs_;                                     ///< Auxiliary vector for assemble methods

    /// Following data members represent Element quantities and FE quantities
    FeQ<Scalar> JxW_;
    FeQ<Scalar> JxW_side_;
    ElQ<Vector> normal_;
    FeQArray<Vector> deform_;
    FeQArray<Vector> deform_side_;
    FeQArray<Tensor> grad_deform_;
    FeQArray<Scalar> div_deform_;
    FeQJoin<Vector> deform_join_;


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
    Rhs_ComputeLocal(EqFields *eq_fields, EqData *eq_data, PatchFEValues<3> *fe_values)
    : Rhs_FullAssembly<dim>(eq_fields, eq_data, fe_values) {}

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
    Rhs_EvalFields(EqFields *eq_fields, EqData *eq_data, PatchFEValues<3> *fe_values)
    : Rhs_FullAssembly<dim>(eq_fields, eq_data, fe_values) {}

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

