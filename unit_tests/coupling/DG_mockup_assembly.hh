#ifndef DG_MOCKUP_ASSEMBLYA_HH_
#define DG_MOCKUP_ASSEMBLYA_HH_

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "DG_mockup.hh"
#include "fem/fe_p.hh"
#include "fem/patch_fe_values.hh"
#include "quadrature/quadrature_lib.hh"
#include "coupling/balance.hh"
#include "fields/field_value_cache.hh"



/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim>
class Mass_FullAssembly : public AssemblyBasePatch<dim>
{
public:
    typedef equation_data::EqFields EqFields;
    typedef equation_data::EqData EqData;

    static constexpr const char * name() { return "MassAssembly"; }

    /// Constructor.
    Mass_FullAssembly(EqFields *eq_fields, EqData *eq_data, PatchFEValues<3> *fe_values)
    : AssemblyBasePatch<dim>(eq_data->dg_order, fe_values), eq_fields_(eq_fields), eq_data_(eq_data),
      JxW_( this->fe_values_->JxW(this->quad_) ),
      conc_shape_( this->fe_values_->scalar_shape(this->quad_) ) {
        this->active_integrals_ = ActiveIntegrals::bulk;
        this->used_fields_ += eq_fields_->mass_matrix_coef;
        this->used_fields_ += eq_fields_->retardation_coef;
    }

    /// Destructor.
    ~Mass_FullAssembly() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        this->element_cache_map_ = element_cache_map;

        UpdateFlags u = update_values | update_JxW_values | update_quadrature_points;
        this->fe_values_->template initialize<dim>(*this->quad_, u);
        ndofs_ = this->fe_values_->n_dofs(dim);
        dof_indices_.resize(ndofs_);
        local_matrix_.resize(4*ndofs_*ndofs_);
        local_retardation_balance_vector_.resize(ndofs_);
        local_mass_balance_vector_.resize(ndofs_);

    }


    /// Assemble integral over element
    inline virtual void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        ASSERT_EQ(cell.dim(), dim).error("Dimension of element mismatch!");

        cell.get_dof_indices(dof_indices_);

        for (unsigned int sbi=0; sbi<eq_data_->n_substances(); ++sbi)
        {
            // assemble the local mass matrix
            for (unsigned int i=0; i<ndofs_; i++)
            {
                for (unsigned int j=0; j<ndofs_; j++)
                {
                    local_matrix_[i*ndofs_+j] = 0;
                    for (auto p : this->bulk_points(element_patch_idx) )
                    {
                        local_matrix_[i*ndofs_+j] += (eq_fields_->mass_matrix_coef(p)+eq_fields_->retardation_coef[sbi](p)) *
                                conc_shape_(j,p)*conc_shape_(i,p)*JxW_(p);
                    }
                }
            }

            for (unsigned int i=0; i<ndofs_; i++)
            {
                local_mass_balance_vector_[i] = 0;
                local_retardation_balance_vector_[i] = 0;
                for (auto p : this->bulk_points(element_patch_idx) )
                {
                    local_mass_balance_vector_[i] += eq_fields_->mass_matrix_coef(p)*conc_shape_(i,p)*JxW_(p);
                    local_retardation_balance_vector_[i] -= eq_fields_->retardation_coef[sbi](p)*conc_shape_(i,p)*JxW_(p);
                }
            }

            this->cell_integral_set_values(sbi);
        }
    }

    protected:
        virtual void cell_integral_set_values(unsigned int sbi) {
            this->eq_data_->ls_dt[sbi]->mat_set_values(this->ndofs_, &(this->dof_indices_[0]), this->ndofs_, &(this->dof_indices_[0]), &(this->local_matrix_[0]));
            VecSetValues(this->eq_data_->ret_vec[sbi], this->ndofs_, &(this->dof_indices_[0]), &(this->local_retardation_balance_vector_[0]), ADD_VALUES);
        }

        /// Data objects shared with DGMockup
        EqFields *eq_fields_;
        EqData *eq_data_;

        /// Sub field set contains fields used in calculation.
        FieldSet used_fields_;

        unsigned int ndofs_;                                      ///< Number of dofs
        vector<LongIdx> dof_indices_;                             ///< Vector of global DOF indices
        vector<PetscScalar> local_matrix_;                        ///< Auxiliary vector for assemble methods
        vector<PetscScalar> local_retardation_balance_vector_;    ///< Auxiliary vector for assemble mass matrix.
        vector<PetscScalar> local_mass_balance_vector_;           ///< Same as previous.

        ElQ<Scalar> JxW_;
        FeQ<Scalar> conc_shape_;

        template < template<IntDim...> class DimAssembly>
        friend class GenericAssembly;

};


template <unsigned int dim>
class Mass_ComputeLocal : public Mass_FullAssembly<dim>
{
public:
    typedef equation_data::EqFields EqFields;
    typedef equation_data::EqData EqData;

    static constexpr const char * name() { return "MassAssembly"; }

    /// Constructor.
    Mass_ComputeLocal(EqFields *eq_fields, EqData *eq_data, PatchFEValues<3> *fe_values)
    : Mass_FullAssembly<dim>(eq_fields, eq_data, fe_values) {}

    /// Destructor.
    ~Mass_ComputeLocal() {}

protected:
    void cell_integral_set_values(unsigned int sbi) override {}

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


template <unsigned int dim>
class Mass_EvalFields : public Mass_FullAssembly<dim>
{
public:
    typedef equation_data::EqFields EqFields;
    typedef equation_data::EqData EqData;

    static constexpr const char * name() { return "MassAssembly"; }

    /// Constructor.
    Mass_EvalFields(EqFields *eq_fields, EqData *eq_data, PatchFEValues<3> *fe_values)
    : Mass_FullAssembly<dim>(eq_fields, eq_data, fe_values) {}

    /// Destructor.
    ~Mass_EvalFields() {}

    void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx) override {}

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim>
class Stiffness_FullAssembly : public AssemblyBasePatch<dim>
{
public:
    typedef equation_data::EqFields EqFields;
    typedef equation_data::EqData EqData;

    static constexpr const char * name() { return "StiffnessAssembly"; }

    /// Constructor.
    Stiffness_FullAssembly(EqFields *eq_fields, EqData *eq_data, PatchFEValues<3> *fe_values)
    : AssemblyBasePatch<dim>(eq_data->dg_order, fe_values), eq_fields_(eq_fields), eq_data_(eq_data),
      JxW_( this->fe_values_->JxW(this->quad_) ),
      JxW_side_( this->fe_values_->JxW_side(this->quad_low_) ),
      normal_( this->fe_values_->normal_vector(this->quad_low_) ),
      conc_shape_( this->fe_values_->scalar_shape(this->quad_) ),
      conc_shape_side_( this->fe_values_->scalar_shape_side(this->quad_low_) ),
      conc_grad_( this->fe_values_->grad_scalar_shape(this->quad_) ),
	  conc_grad_side_( this->fe_values_->grad_scalar_shape_side(this->quad_low_) ),
      conc_join_shape_(make_iter<JoinShapeAccessor<Scalar>>(JoinShapeAccessor<Scalar>()), make_iter<JoinShapeAccessor<Scalar>>(JoinShapeAccessor<Scalar>())) {
        this->active_integrals_ = (ActiveIntegrals::bulk | ActiveIntegrals::edge | ActiveIntegrals::coupling | ActiveIntegrals::boundary);
        this->used_fields_ += eq_fields_->advection_coef;
        this->used_fields_ += eq_fields_->diffusion_coef;
        this->used_fields_ += eq_fields_->cross_section;
        this->used_fields_ += eq_fields_->dg_penalty;
        this->used_fields_ += eq_fields_->sources_sigma_out;
        this->used_fields_ += eq_fields_->fracture_sigma;
        this->used_fields_ += eq_fields_->bc_type;
        this->used_fields_ += eq_fields_->bc_robin_sigma;
    }

    /// Destructor.
    ~Stiffness_FullAssembly() {
        for (auto a : averages) if (a != nullptr) delete[] a;
        for (auto a : waverages) if (a != nullptr) delete[] a;
        for (auto a : jumps) if (a != nullptr) delete[] a;
    }

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        this->element_cache_map_ = element_cache_map;

        UpdateFlags u = update_values | update_gradients | update_JxW_values | update_quadrature_points;
        UpdateFlags u_side = update_values | update_gradients | update_side_JxW_values | update_normal_vectors | update_quadrature_points;
        this->fe_values_->template initialize<dim>(*this->quad_, u);
        this->fe_values_->template initialize<dim>(*this->quad_low_, u_side);
        if (dim>1)
            conc_join_shape_ = Range< JoinShapeAccessor<Scalar> >( this->fe_values_->scalar_join_shape( {this->quad_low_, this->quad_low_} ) );
        if (dim==1) { // print to log only one time
            DebugOut() << "List of StiffnessAssemblyDG FEValues (cell) updates flags: " << this->print_update_flags(u);
            DebugOut() << "List of StiffnessAssemblyDG FEValues (side) updates flags: " << this->print_update_flags(u_side);
        }
        ndofs_ = this->fe_values_->n_dofs(dim);
        qsize_lower_dim_ = this->quad_low_->size();
        dof_indices_.resize(ndofs_);
        side_dof_indices_vb_.resize(2*ndofs_);
        local_matrix_.resize(4*ndofs_*ndofs_);

        for (unsigned int sid=0; sid<eq_data_->max_edg_sides; sid++)
        {
            side_dof_indices_.push_back( vector<LongIdx>(ndofs_) );
        }

        averages.resize(eq_data_->max_edg_sides);
        for (unsigned int s=0; s<eq_data_->max_edg_sides; s++)
            averages[s] = new double[qsize_lower_dim_*ndofs_];
        waverages.resize(2);
        jumps.resize(2);
        for (unsigned int s=0; s<2; s++)
        {
            waverages[s] = new double[qsize_lower_dim_*ndofs_];
            jumps[s] = new double[qsize_lower_dim_*ndofs_];
        }
    }


    /// Assembles the cell (volume) integral into the stiffness matrix.
    inline virtual void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        ASSERT_EQ(cell.dim(), dim).error("Dimension of element mismatch!");
        if (!cell.is_own()) return;

        cell.get_dof_indices(dof_indices_);

        // assemble the local stiffness matrix
        for (unsigned int sbi=0; sbi<eq_data_->n_substances(); sbi++)
        {
            for (unsigned int i=0; i<ndofs_; i++)
                for (unsigned int j=0; j<ndofs_; j++)
                    local_matrix_[i*ndofs_+j] = 0;

            for (auto p : this->bulk_points(element_patch_idx) )
            {
                for (unsigned int i=0; i<ndofs_; i++)
                {
                    arma::vec3 Kt_grad_i = eq_fields_->diffusion_coef[sbi](p).t()*conc_grad_(i,p);
                    double ad_dot_grad_i = arma::dot(eq_fields_->advection_coef[sbi](p), conc_grad_(i,p));

                    for (unsigned int j=0; j<ndofs_; j++)
                        local_matrix_[i*ndofs_+j] += (arma::dot(Kt_grad_i, conc_grad_(j,p))
                                                  -conc_shape_(j,p)*ad_dot_grad_i
                                                  +eq_fields_->sources_sigma_out[sbi](p)*conc_shape_(j,p)*conc_shape_(i,p))*JxW_(p);
                }
            }
            this->cell_integral_set_values(sbi);
        }
    }


    /// Assembles the fluxes on the boundary.
    inline virtual void boundary_side_integral(DHCellSide cell_side)
    {
        ASSERT_EQ(cell_side.dim(), dim).error("Dimension of element mismatch!");
        if (!cell_side.cell().is_own()) return;

        Side side = cell_side.side();
        const DHCellAccessor &cell = cell_side.cell();

        cell.get_dof_indices(dof_indices_);
        unsigned int k;
        double gamma_l;

        for (unsigned int sbi=0; sbi<eq_data_->n_substances(); sbi++)
        {
            std::fill(local_matrix_.begin(), local_matrix_.end(), 0);

            // On Neumann boundaries we have only term from integrating by parts the advective term,
            // on Dirichlet boundaries we additionally apply the penalty which enforces the prescribed value.
            double side_flux = 0;
            k=0;
            for (auto p : this->boundary_points(cell_side) ) {
                side_flux += arma::dot(eq_fields_->advection_coef[sbi](p), normal_(p))*JxW_side_(p);
                k++;
            }
            double transport_flux = side_flux/side.measure();

            auto p_side = *( this->boundary_points(cell_side).begin() );
            auto p_bdr = p_side.point_bdr( side.cond().element_accessor() );
            unsigned int bc_type = eq_fields_->bc_type[sbi](p_bdr);
            if (bc_type == DGMockup<Mass_FullAssembly, Stiffness_FullAssembly, Sources_FullAssembly>::abc_dirichlet)
            {
                // set up the parameters for DG method
                k=0; // temporary solution, set dif_coef until set_DG_parameters_boundary will not be removed
                for (auto p : this->boundary_points(cell_side) ) {
                    eq_data_->dif_coef[sbi][k] = eq_fields_->diffusion_coef[sbi](p);
                    k++;
                }
                auto p = *( this->boundary_points(cell_side).begin() );
                eq_data_->set_DG_parameters_boundary(side, qsize_lower_dim_, eq_data_->dif_coef[sbi], transport_flux, normal_(p), eq_fields_->dg_penalty[sbi](p_side), gamma_l);
                transport_flux += gamma_l;
            }

            // fluxes and penalty
            k=0;
            for (auto p : this->boundary_points(cell_side) )
            {
                double flux_times_JxW;
                if (bc_type == DGMockup<Mass_FullAssembly, Stiffness_FullAssembly, Sources_FullAssembly>::abc_total_flux)
                {
                    //sigma_ corresponds to robin_sigma
                    auto p_bdr = p.point_bdr(side.cond().element_accessor());
                    flux_times_JxW = eq_fields_->cross_section(p)*eq_fields_->bc_robin_sigma[sbi](p_bdr)*JxW_side_(p);
                }
                else if (bc_type == DGMockup<Mass_FullAssembly, Stiffness_FullAssembly, Sources_FullAssembly>::abc_diffusive_flux)
                {
                    auto p_bdr = p.point_bdr(side.cond().element_accessor());
                    flux_times_JxW = (transport_flux + eq_fields_->cross_section(p)*eq_fields_->bc_robin_sigma[sbi](p_bdr))*JxW_side_(p);
                }
                else if (bc_type == DGMockup<Mass_FullAssembly, Stiffness_FullAssembly, Sources_FullAssembly>::abc_inflow && side_flux < 0)
                    flux_times_JxW = 0;
                else
                    flux_times_JxW = transport_flux*JxW_side_(p);

                for (unsigned int i=0; i<ndofs_; i++)
                {
                    for (unsigned int j=0; j<ndofs_; j++)
                    {
                        // flux due to advection and penalty
                        local_matrix_[i*ndofs_+j] += flux_times_JxW*conc_shape_side_(i,p)*conc_shape_side_(j,p);

                        // flux due to diffusion (only on dirichlet and inflow boundary)
                        if (bc_type == DGMockup<Mass_FullAssembly, Stiffness_FullAssembly, Sources_FullAssembly>::abc_dirichlet)
                            local_matrix_[i*ndofs_+j] -= (arma::dot(eq_fields_->diffusion_coef[sbi](p)*conc_grad_side_(j,p),normal_(p))*conc_shape_side_(i,p)
                                    + arma::dot(eq_fields_->diffusion_coef[sbi](p)*conc_grad_side_(i,p),normal_(p))*conc_shape_side_(j,p)*eq_data_->dg_variant
                                    )*JxW_side_(p);
                    }
                }
                k++;
            }

            this->boundary_side_integral_set_values(sbi);
        }
    }


    /// Assembles the fluxes between sides of elements of the same dimension.
    inline virtual void edge_integral(RangeConvert<DHEdgeSide, DHCellSide> edge_side_range) {
        ASSERT_EQ(edge_side_range.begin()->element().dim(), dim).error("Dimension of element mismatch!");

        unsigned int k;
        double gamma_l, omega[2], transport_flux, delta[2], delta_sum;
        double aniso1, aniso2;
        double local_alpha=0.0;
        int sid=0, s1, s2;
        for( DHCellSide edge_side : edge_side_range )
        {
            auto dh_edge_cell = eq_data_->dh_->cell_accessor_from_element( edge_side.elem_idx() );
            dh_edge_cell.get_dof_indices(side_dof_indices_[sid]);
            ++sid;
        }
        auto zero_edge_side = *edge_side_range.begin();
        auto p = *( this->edge_points(zero_edge_side).begin() );
        arma::vec3 normal_vector = normal_(p);

        // fluxes and penalty
        for (unsigned int sbi=0; sbi<eq_data_->n_substances(); sbi++)
        {
            vector<double> fluxes(edge_side_range.begin()->n_edge_sides());
            double pflux = 0, nflux = 0; // calculate the total in- and out-flux through the edge
            sid=0;
            for( DHCellSide edge_side : edge_side_range )
            {
                fluxes[sid] = 0;
                for (auto p : this->edge_points(edge_side) ) {
                    fluxes[sid] += arma::dot(eq_fields_->advection_coef[sbi](p), normal_(p))*JxW_side_(p);
                }
                fluxes[sid] /= edge_side.measure();
                if (fluxes[sid] > 0)
                    pflux += fluxes[sid];
                else
                    nflux += fluxes[sid];
                ++sid;
            }

            // precompute averages of shape functions over pairs of sides
            s1=0;
            for (DHCellSide edge_side : edge_side_range)
            {
                k=0;
                for (auto p : this->edge_points(edge_side) )
                {
                    for (unsigned int i=0; i<ndofs_; i++)
                        averages[s1][k*ndofs_+i] = conc_shape_side_(i,p)*0.5;
                    k++;
                }
                s1++;
            }



            s1=0;
            for( DHCellSide edge_side1 : edge_side_range )
            {
                s2=-1; // need increment at begin of loop (see conditionally 'continue' directions)
                for( DHCellSide edge_side2 : edge_side_range )
                {
                    s2++;
                    if (s2<=s1) continue;
                    ASSERT(edge_side1.is_valid()).error("Invalid side of edge.");

                    auto p = *( this->edge_points(edge_side1).begin() );
                    arma::vec3 nv = normal_(p);

                    // set up the parameters for DG method
                    // calculate the flux from edge_side1 to edge_side2
                    if (fluxes[s2] > 0 && fluxes[s1] < 0)
                        transport_flux = fluxes[s1]*fabs(fluxes[s2]/pflux);
                    else if (fluxes[s2] < 0 && fluxes[s1] > 0)
                        transport_flux = fluxes[s1]*fabs(fluxes[s2]/nflux);
                    else
                        transport_flux = 0;

                    gamma_l = 0.5*fabs(transport_flux);

                    delta[0] = 0;
                    delta[1] = 0;
                    for (auto p1 : this->edge_points(edge_side1) )
                    {
                        auto p2 = p1.point_on(edge_side2);
                        delta[0] += dot(eq_fields_->diffusion_coef[sbi](p1)*normal_vector,normal_vector);
                        delta[1] += dot(eq_fields_->diffusion_coef[sbi](p2)*normal_vector,normal_vector);
                        local_alpha = max(eq_fields_->dg_penalty[sbi](p1), eq_fields_->dg_penalty[sbi](p2));
                    }
                    delta[0] /= qsize_lower_dim_;
                    delta[1] /= qsize_lower_dim_;

                    delta_sum = delta[0] + delta[1];

//                        if (delta_sum > numeric_limits<double>::epsilon())
                    if (fabs(delta_sum) > 0)
                    {
                        omega[0] = delta[1]/delta_sum;
                        omega[1] = delta[0]/delta_sum;
                        double h = edge_side1.diameter();
                        aniso1 = eq_data_->elem_anisotropy(edge_side1.element());
                        aniso2 = eq_data_->elem_anisotropy(edge_side2.element());
                        gamma_l += local_alpha/h*aniso1*aniso2*(delta[0]*delta[1]/delta_sum);
                    }
                    else
                        for (int i=0; i<2; i++) omega[i] = 0;
                    // end of set up the parameters for DG method

                    int sd[2]; bool is_side_own[2];
                    sd[0] = s1; is_side_own[0] = edge_side1.cell().is_own();
                    sd[1] = s2; is_side_own[1] = edge_side2.cell().is_own();

                    // precompute jumps and weighted averages of shape functions over the pair of sides (s1,s2)
                    k=0;
                    for (auto p1 : this->edge_points(edge_side1) )
                    {
                        auto p2 = p1.point_on(edge_side2);
                        for (unsigned int i=0; i<ndofs_; i++)
                        {
                            jumps[0][k*ndofs_+i] = conc_shape_side_(i,p1);
                            jumps[1][k*ndofs_+i] = - conc_shape_side_(i,p2);
                            waverages[0][k*ndofs_+i] = arma::dot(eq_fields_->diffusion_coef[sbi](p1)*conc_grad_side_(i,p1),nv)*omega[0];
                            waverages[1][k*ndofs_+i] = arma::dot(eq_fields_->diffusion_coef[sbi](p2)*conc_grad_side_(i,p2),nv)*omega[1];
                        }
                        k++;
                    }

                    // For selected pair of elements:
                    for (int n=0; n<2; n++)
                    {
                        if (!is_side_own[n]) continue;

                        for (int m=0; m<2; m++)
                        {
                            for (unsigned int i=0; i<this->fe_values_->n_dofs(dim); i++)
                                for (unsigned int j=0; j<this->fe_values_->n_dofs(dim); j++)
                                    local_matrix_[i*this->fe_values_->n_dofs(dim)+j] = 0;

                            k=0;
                            for (auto p1 : this->edge_points(zero_edge_side) )
                            //for (k=0; k<this->quad_low_->size(); ++k)
                            {
                                for (unsigned int i=0; i<this->fe_values_->n_dofs(dim); i++)
                                {
                                    for (unsigned int j=0; j<this->fe_values_->n_dofs(dim); j++)
                                    {
                                        int index = i*this->fe_values_->n_dofs(dim)+j;

                                        local_matrix_[index] += (
                                            // flux due to transport (applied on interior edges) (average times jump)
                                            transport_flux*jumps[n][k*ndofs_+i]*averages[sd[m]][k*ndofs_+j]

                                            // penalty enforcing continuity across edges (applied on interior and Dirichlet edges) (jump times jump)
                                            + gamma_l*jumps[n][k*ndofs_+i]*jumps[m][k*ndofs_+j]

                                        // terms due to diffusion
                                            - jumps[n][k*ndofs_+i]*waverages[m][k*ndofs_+j]
                                            - eq_data_->dg_variant*waverages[n][k*ndofs_+i]*jumps[m][k*ndofs_+j]
                                            )*JxW_side_(p1) + LocalSystem::almost_zero;
                                    }
                                }
                                k++;
                            }
                            this->edge_integral_set_values(sbi, ndofs_, m, n, sd);
                        }
                    }
                }
            s1++;
            }
        }
    }


    /// Assembles the fluxes between elements of different dimensions.
    inline virtual void dimjoin_intergral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) {
        if (dim == 1) return;
        ASSERT_EQ(cell_lower_dim.dim(), dim-1).error("Dimension of element mismatch!");

        // Note: use data members csection_ and velocity_ for appropriate quantities of lower dim element

        unsigned int n_dofs[2];
        unsigned int n_indices = cell_lower_dim.get_dof_indices(dof_indices_);
        for(unsigned int i=0; i<n_indices; ++i) {
            side_dof_indices_vb_[i] = dof_indices_[i];
        }
        n_dofs[0] = conc_join_shape_.begin()->n_dofs_low();

        DHCellAccessor cell_higher_dim = eq_data_->dh_->cell_accessor_from_element( neighb_side.element().idx() );
        n_indices = cell_higher_dim.get_dof_indices(dof_indices_);
        for(unsigned int i=0; i<n_indices; ++i) {
            side_dof_indices_vb_[i+n_dofs[0]] = dof_indices_[i];
        }
        n_dofs[1] = conc_join_shape_.begin()->n_dofs_high();

        // Testing element if they belong to local partition.
        bool own_element_id[2];
        own_element_id[0] = cell_lower_dim.is_own();
        own_element_id[1] = cell_higher_dim.is_own();

        for (unsigned int sbi=0; sbi<eq_data_->n_substances(); sbi++) // Optimize: SWAP LOOPS
        {
            for (unsigned int i=0; i<n_dofs[0]+n_dofs[1]; i++)
                for (unsigned int j=0; j<n_dofs[0]+n_dofs[1]; j++)
                    local_matrix_[i*(n_dofs[0]+n_dofs[1])+j] = 0;

            // set transmission conditions
            for (auto p_high : this->coupling_points(neighb_side) )
            {
                auto p_low = p_high.lower_dim(cell_lower_dim);
                auto p_low_conv = p_high.lower_dim_converted(cell_lower_dim);
                // The communication flux has two parts:
                // - "diffusive" term containing sigma
                // - "advective" term representing usual upwind
                //
                // The calculation differs from the reference manual, since ad_coef and dif_coef have different meaning
                // than b and A in the manual.
                // In calculation of sigma there appears one more csection_lower in the denominator.
                double sigma = eq_fields_->fracture_sigma[sbi](p_low)*arma::dot(eq_fields_->diffusion_coef[sbi](p_low)*normal_(p_high),normal_(p_high))*
                        2*eq_fields_->cross_section(p_high)*eq_fields_->cross_section(p_high)/(eq_fields_->cross_section(p_low)*eq_fields_->cross_section(p_low));

                double transport_flux = arma::dot(eq_fields_->advection_coef[sbi](p_high), normal_(p_high));

                for( auto conc_shape_i : conc_join_shape_) {
                    uint is_high_i = conc_shape_i.is_high_dim();
                    if (!own_element_id[is_high_i]) continue;
                    uint i_mat_idx = conc_shape_i.join_idx(); // i + is_high * n_dofs_low
                    double diff_shape_i = conc_shape_i(p_high) - conc_shape_i(p_low_conv);
                    for( auto conc_shape_j : conc_join_shape_) {
                        uint j_mat_idx = conc_shape_j.join_idx();
                        local_matrix_[i_mat_idx * (n_dofs[0]+n_dofs[1]) + j_mat_idx] += (
                                sigma * diff_shape_i * (conc_shape_j(p_high) - conc_shape_j(p_low_conv))
                                + diff_shape_i * ( max(0.,transport_flux) * conc_shape_j(p_high) + min(0.,transport_flux) * conc_shape_j(p_low_conv))
						    )*JxW_side_(p_high) + LocalSystem::almost_zero;
                    }
                }
            }
            this->dimjoin_intergral_set_values(sbi, n_dofs);
        }
    }


protected:
    /// Holds set of element_patch_idx and side_idx of processed edge.
    struct EdgeValuesMap {
        std::vector<unsigned int> element_patch_idx_;      ///< Index of element in patch
        std::vector<unsigned int> side_idx_;               ///< Index of side in element
    };

    virtual void cell_integral_set_values(unsigned int sbi) {
        this->eq_data_->ls[sbi]->mat_set_values(this->ndofs_, &(this->dof_indices_[0]), this->ndofs_, &(this->dof_indices_[0]), &(this->local_matrix_[0]));
    }

    virtual void boundary_side_integral_set_values(unsigned int sbi) {
        this->eq_data_->ls[sbi]->mat_set_values(this->ndofs_, &(this->dof_indices_[0]), this->ndofs_, &(this->dof_indices_[0]), &(this->local_matrix_[0]));
    }

    virtual void edge_integral_set_values(unsigned int sbi, unsigned int n_dofs, int m, int n, int sd[2]) {
        this->eq_data_->ls[sbi]->mat_set_values(n_dofs, &(this->side_dof_indices_[sd[n]][0]), n_dofs, &(this->side_dof_indices_[sd[m]][0]), &(this->local_matrix_[0]));
    }

    virtual void dimjoin_intergral_set_values(unsigned int sbi, unsigned int n_dofs[2]) {
        this-> eq_data_->ls[sbi]->mat_set_values(n_dofs[0]+n_dofs[1], &(this->side_dof_indices_vb_[0]), n_dofs[0]+n_dofs[1], &(this->side_dof_indices_vb_[0]), &(this->local_matrix_[0]));
    }


    /// Data objects shared with TransportDG
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int ndofs_;                                      ///< Number of dofs
    unsigned int qsize_lower_dim_;                            ///< Size of quadrature of dim-1

    vector<LongIdx> dof_indices_;                             ///< Vector of global DOF indices
    vector< vector<LongIdx> > side_dof_indices_;              ///< Vector of vectors of side DOF indices
    vector<LongIdx> side_dof_indices_vb_;                     ///< Vector of side DOF indices (assemble element-side fluxex)
    vector<PetscScalar> local_matrix_;                        ///< Auxiliary vector for assemble methods

    vector<double*> averages;                                 ///< Auxiliary storage for averages of shape functions.
    vector<double*> waverages;                                ///< Auxiliary storage for weighted averages of shape functions.
    vector<double*> jumps;                                    ///< Auxiliary storage for jumps of shape functions.

    ElQ<Scalar> JxW_;
    ElQ<Scalar> JxW_side_;
    ElQ<Vector> normal_;
    FeQ<Scalar> conc_shape_;
    FeQ<Scalar> conc_shape_side_;
    FeQ<Vector> conc_grad_;
    FeQ<Vector> conc_grad_side_;
    Range< JoinShapeAccessor<Scalar> > conc_join_shape_;

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
    void cell_integral_set_values(unsigned int sbi) override {}

    void boundary_side_integral_set_values(unsigned int sbi) override {}

    void edge_integral_set_values(unsigned int sbi, unsigned int n_dofs, int m, int n, int sd[2]) override {}

    void dimjoin_intergral_set_values(unsigned int sbi, unsigned int n_dofs[2]) override {}

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


    /// Assembles the fluxes between sides of elements of the same dimension.
    inline void edge_integral(RangeConvert<DHEdgeSide, DHCellSide> edge_side_range) override {}


    /// Assembles the fluxes between elements of different dimensions.
    inline void dimjoin_intergral(DHCellAccessor cell_lower_dim, DHCellSide neighb_side) override {}

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


/**
 * Auxiliary container class for Finite element and related objects of given dimension.
 */
template <unsigned int dim>
class Sources_FullAssembly : public AssemblyBasePatch<dim>
{
public:
    typedef equation_data::EqFields EqFields;
    typedef equation_data::EqData EqData;

    static constexpr const char * name() { return "SourcesAssembly"; }

    /// Constructor.
    Sources_FullAssembly(EqFields *eq_fields, EqData *eq_data, PatchFEValues<3> *fe_values)
    : AssemblyBasePatch<dim>(eq_data->dg_order, fe_values), eq_fields_(eq_fields), eq_data_(eq_data),
      JxW_( this->fe_values_->JxW(this->quad_) ),
      conc_shape_( this->fe_values_->scalar_shape(this->quad_) ) {
        this->active_integrals_ = ActiveIntegrals::bulk;
        this->used_fields_ += eq_fields_->sources_density_out;
        this->used_fields_ += eq_fields_->sources_conc_out;
        this->used_fields_ += eq_fields_->sources_sigma_out;
    }

    /// Destructor.
    ~Sources_FullAssembly() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        this->element_cache_map_ = element_cache_map;

        UpdateFlags u = update_values | update_JxW_values | update_quadrature_points;
        this->fe_values_->template initialize<dim>(*this->quad_, u);
        ndofs_ = this->fe_values_->n_dofs(dim);
        dof_indices_.resize(ndofs_);
        local_rhs_.resize(ndofs_);
        local_source_balance_vector_.resize(ndofs_);
        local_source_balance_rhs_.resize(ndofs_);
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        ASSERT_EQ(cell.dim(), dim).error("Dimension of element mismatch!");

        ElementAccessor<3> elm = cell.elm();
        double source;

        cell.get_dof_indices(dof_indices_);

        // assemble the local stiffness matrix
        for (unsigned int sbi=0; sbi<eq_data_->n_substances(); sbi++)
        {
            fill_n( &(local_rhs_[0]), ndofs_, 0 );
            local_source_balance_vector_.assign(ndofs_, 0);
            local_source_balance_rhs_.assign(ndofs_, 0);

            for (auto p : this->bulk_points(element_patch_idx) )
            {
                source = (eq_fields_->sources_density_out[sbi](p) + eq_fields_->sources_conc_out[sbi](p)*eq_fields_->sources_sigma_out[sbi](p))*JxW_(p);

                for (unsigned int i=0; i<ndofs_; i++)
                    local_rhs_[i] += source*conc_shape_(i,p);
            }
            this->cell_integral_set_values(sbi);

            for (unsigned int i=0; i<ndofs_; i++)
            {
                for (auto p : this->bulk_points(element_patch_idx) )
                {
                    local_source_balance_vector_[i] -= eq_fields_->sources_sigma_out[sbi](p)*conc_shape_(i,p)*JxW_(p);
                }

                local_source_balance_rhs_[i] += local_rhs_[i];
            }
        }
    }


protected:
    virtual void cell_integral_set_values(unsigned int sbi) {
        this->eq_data_->ls[sbi]->rhs_set_values(this->ndofs_, &(this->dof_indices_[0]), &(this->local_rhs_[0]));
    }


    /// Data objects shared with TransportDG
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    unsigned int ndofs_;                                      ///< Number of dofs

    vector<LongIdx> dof_indices_;                             ///< Vector of global DOF indices
    vector<PetscScalar> local_rhs_;                           ///< Auxiliary vector for set_sources method.
    vector<PetscScalar> local_source_balance_vector_;         ///< Auxiliary vector for set_sources method.
    vector<PetscScalar> local_source_balance_rhs_;            ///< Auxiliary vector for set_sources method.

    ElQ<Scalar> JxW_;
    FeQ<Scalar> conc_shape_;

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


template <unsigned int dim>
class Sources_ComputeLocal : public Sources_FullAssembly<dim>
{
public:
    typedef equation_data::EqFields EqFields;
    typedef equation_data::EqData EqData;

    static constexpr const char * name() { return "SourcesAssembly"; }

    /// Constructor.
    Sources_ComputeLocal(EqFields *eq_fields, EqData *eq_data, PatchFEValues<3> *fe_values)
    : Sources_FullAssembly<dim>(eq_fields, eq_data, fe_values) {}

    /// Destructor.
    ~Sources_ComputeLocal() {}

protected:
    virtual void cell_integral_set_values(unsigned int sbi) {}

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


template <unsigned int dim>
class Sources_EvalFields : public Sources_FullAssembly<dim>
{
public:
    typedef equation_data::EqFields EqFields;
    typedef equation_data::EqData EqData;

    static constexpr const char * name() { return "SourcesAssembly"; }

    /// Constructor.
    Sources_EvalFields(EqFields *eq_fields, EqData *eq_data, PatchFEValues<3> *fe_values)
    : Sources_FullAssembly<dim>(eq_fields, eq_data, fe_values) {}

    /// Destructor.
    ~Sources_EvalFields() {}

    void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx) override {}

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


#endif /* DG_MOCKUP_ASSEMBLYA_HH_ */

