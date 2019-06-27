/*
 * assembly_lmh.hh
 *
 *  Created on: Apr 30, 2016
 *      Author: jb
 */

#ifndef SRC_FLOW_ASSEMBLY_LMH_HH_
#define SRC_FLOW_ASSEMBLY_LMH_HH_

#include "mesh/long_idx.hh"
#include "flow/darcy_flow_assembly.hh"
#include "soil_models.hh"
#include "coupling/balance.hh"

#include "badiff.h"


/**
 * Prove of concept for general assembly scheme.
 * Ideas:
 * - Local assembly class for given DIM provides assembly methods for:
 *   - local assembly on elements of dimension DIM
 *   - local assembly on edges/faces of neighbouring elemnts of dimension DIM
 *   - local assembly of interactiong elements DIM and OTHER_DIM with OTHER_DIM < DIM
 * - Such class template would be necessary to create particular instance of Assembler class, that
 *   takes job of pass through the mesh, precomputing necessary fields, etc.
 *   This allows to manage order of assembled entities in arbitrary way in order to use caches and vectorisation efficiently.
 * - Separate local assembly template would be necessary for every pass through the mesh, but it may delegate its actions
 *   to a single general local assembly class.
 * - The local assembly class gets an Accessor object for particular domain of the assembly, i.e. ElementAccessor for element assembly,
 *   FaceAccessor for interface intergrals assembly, and InteractionAccessor for assembly of interaction of two elements, e.g. for XFEM
 *   there may be an interaction between every 3D element within the enrichment with appropriate 1D elemnt on the well.
 *   These accessors provides methods to access fields values as well as the DOF indices on the element(s).
 * - Results of the assembly are passed to the global linear algebra objects collected in the RichardSystem class, global indices
 *   (which are still local indicies of the single process) provided by accessors are used.
 *
 *
 * TODO:
 * - finish functional Richards model
 * - move whole internals of the assembly_mh_matrix into local assembly classes
 * - mean while: imporve accessors, imporve mesh and DOF handler interface; possibly create new mesh implementation used in Darcy first and then
 *   apply it to other equations
 */


template<int dim>
class AssemblyLMH : public AssemblyMH<dim> {
public:

    typedef std::shared_ptr<RichardsLMH::EqData> AssemblyDataPtr;

    AssemblyLMH(AssemblyDataPtr data)
    : AssemblyMH<dim>(data),
      ad_(data),
      genuchten_on(false),
      cross_section(1.0),
      soil_model(data->soil_model_),
      cr_disc_dofs(dim+1)
    {}

protected:
    void reset_soil_model(const DHCellAccessor& dh_cell) {
        ElementAccessor<3> ele = dh_cell.elm();
        genuchten_on = (this->ad_->genuchten_p_head_scale.field_result({ele.region()}) != result_zeros);
        if (genuchten_on) {
            SoilData soil_data;
            soil_data.n =  this->ad_->genuchten_n_exponent.value(ele.centre(), ele);
            soil_data.alpha = this->ad_->genuchten_p_head_scale.value(ele.centre(), ele);
            soil_data.Qr = this->ad_->water_content_residual.value(ele.centre(), ele);
            soil_data.Qs = this->ad_->water_content_saturated.value(ele.centre(), ele);
            soil_data.Ks = this->ad_->conductivity.value(ele.centre(), ele);
            //soil_data.cut_fraction = 0.99; // set by model

            soil_model->reset(soil_data);
        }
    }

    void update_water_content(const DHCellAccessor& dh_cell) {
        reset_soil_model(dh_cell);
        ElementAccessor<3> ele = dh_cell.elm();
        double storativity = this->ad_->storativity.value(ele.centre(), ele);

        for (unsigned int i=0; i<ele->n_sides(); i++) {
            double capacity = 0;
            double water_content = 0;
            double phead = ad_->phead_edge_[ this->edge_indices_[i] ];
            if (genuchten_on) {

                  fadbad::B<double> x_phead(phead);
                  fadbad::B<double> evaluated( soil_model->water_content_diff(x_phead) );
                  evaluated.diff(0,1);
                  water_content = evaluated.val();
                  capacity = x_phead.d(0);
            }
            ad_->capacity[ cr_disc_dofs[i] ] = capacity + storativity;
            ad_->water_content_previous_it[ cr_disc_dofs[i] ] = water_content + storativity * phead;
        }
    }

    void assemble_sides(LocalElementAccessorBase<3> ele) override
    {
        reset_soil_model(ele.dh_cell());
        cross_section = this->ad_->cross_section.value(ele.centre(), ele.element_accessor());


        double conductivity, head;
        if (genuchten_on) {
            conductivity=0;
            head=0;
            ele.dh_cell().cell_with_other_dh(ad_->dh_cr_.get()).get_loc_dof_indices(this->edge_indices_);
            for (unsigned int i=0; i<ele.element_accessor()->n_sides(); i++)
            {
                double phead = ad_->phead_edge_[ this->edge_indices_[i] ];
                conductivity += soil_model->conductivity(phead);
                head += ad_->phead_edge_[ this->edge_indices_[i] ];
            }
            conductivity /= ele.n_sides();
            head /= ele.n_sides();
        } else {
            conductivity = this->ad_->conductivity.value(ele.centre(), ele.element_accessor());
        }

        double scale = 1 / cross_section / conductivity;
        this->assemble_sides_scale(ele,scale);
    }

    /***
     * Called from assembly_local_matrix, assumes precomputed:
     * cross_section, genuchten_on, soil_model
     */
    void assemble_source_term(LocalElementAccessorBase<3> ele) override
    {
        update_dofs(ele.dh_cell());
        update_water_content(ele.dh_cell());
        
        // set lumped source
        double diagonal_coef = ele.measure() * cross_section / ele.n_sides();
        double source_diagonal = diagonal_coef * this->ad_->water_source_density.value(ele.centre(), ele.element_accessor());

        for (unsigned int i=0; i<ele.element_accessor()->n_sides(); i++)
        {

            uint local_side = cr_disc_dofs[i];
            
            double capacity = this->ad_->capacity[local_side];
            double water_content_diff = -ad_->water_content_previous_it[local_side] + ad_->water_content_previous_time[local_side];
            double mass_diagonal = diagonal_coef * capacity;

            /*
            DebugOut().fmt("w diff: {:g}  mass: {:g} w prev: {:g} w time: {:g} c: {:g} p: {:g} z: {:g}",
                    water_content_diff,
                    mass_diagonal * ad_->phead_edge_[local_edge],
                    -ad_->water_content_previous_it[local_side],
                    ad_->water_content_previous_time[local_side],
                    capacity,
                    ad_->phead_edge_[local_edge],
                    ele.centre()[0] );
            */


                double mass_rhs = mass_diagonal * ad_->phead_edge_[ this->edge_indices_[i] ] / this->ad_->time_step_
                                + diagonal_coef * water_content_diff / this->ad_->time_step_;

//                 DBGCOUT(<< "source [" << loc_system_.row_dofs[this->loc_edge_dofs[i]] << ", " << loc_system_.row_dofs[this->loc_edge_dofs[i]]
//                             << "] mat: " << -mass_diagonal/this->ad_->time_step_
//                             << " rhs: " << -source_diagonal - mass_rhs
//                             << "\n");
            this->loc_system_.add_value(this->loc_edge_dofs[i], this->loc_edge_dofs[i],
                                        -mass_diagonal/this->ad_->time_step_,
                                        -source_diagonal - mass_rhs);

            if (ad_->balance != nullptr) {
                ad_->balance->add_mass_vec_value(ad_->water_balance_idx, ele.region().bulk_idx(),
                        diagonal_coef*ad_->water_content_previous_it[local_side]);
                ad_->balance->add_source_values(ad_->water_balance_idx, ele.region().bulk_idx(), {(LongIdx)this->edge_indices_[i]},{0},{source_diagonal});
            }
        }

    }

    /// Updates DoFs for edge pressure vector (dh CR) and for water content vector (dh CR_disc)
    /// Be sure to call it before @p update_water_content().
    void update_dofs(const DHCellAccessor& dh_cell)
    {
        dh_cell.cell_with_other_dh(ad_->dh_cr_.get()).get_loc_dof_indices(this->edge_indices_);
        dh_cell.cell_with_other_dh(ad_->dh_cr_disc_.get()).get_loc_dof_indices(cr_disc_dofs);
    }
    
    void postprocess_velocity_specific(const DHCellAccessor& dh_cell, double edge_scale, double edge_source_term) override
    {   
        update_dofs(dh_cell);
        update_water_content(dh_cell);
        
        dh_cell.get_loc_dof_indices(this->indices_);
        ElementAccessor<3> ele = dh_cell.elm();
        
        for (unsigned int i=0; i<ele->n_sides(); i++) {
            
            double water_content = ad_->water_content_previous_it[ cr_disc_dofs[i] ];
            double water_content_previous_time = ad_->water_content_previous_time[ cr_disc_dofs[i] ];
            
            ad_->data_vec_[this->indices_[this->loc_side_dofs[i]]]
                += edge_source_term - edge_scale * (water_content - water_content_previous_time) / this->ad_->time_step_;
        }
    }

    AssemblyDataPtr ad_;

    bool genuchten_on;
    double cross_section;
    std::shared_ptr<SoilModelBase> soil_model;
    std::vector<LongIdx> cr_disc_dofs;
};



#endif /* SRC_FLOW_ASSEMBLY_LMH_HH_ */
