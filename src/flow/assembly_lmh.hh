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

#include "tools/include_fadbad.hh"  // for "badiff.h"


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
      soil_model(data->soil_model_)
    {}

    void reset_soil_model(LocalElementAccessorBase<3> ele) {
        genuchten_on = (this->ad_->genuchten_p_head_scale.field_result({ele.region()}) != result_zeros);
        if (genuchten_on) {
            SoilData soil_data;
            soil_data.n =  this->ad_->genuchten_n_exponent.value(ele.centre(), ele.element_accessor());
            soil_data.alpha = this->ad_->genuchten_p_head_scale.value(ele.centre(), ele.element_accessor());
            soil_data.Qr = this->ad_->water_content_residual.value(ele.centre(), ele.element_accessor());
            soil_data.Qs = this->ad_->water_content_saturated.value(ele.centre(), ele.element_accessor());
            soil_data.Ks = this->ad_->conductivity.value(ele.centre(), ele.element_accessor());
            //soil_data.cut_fraction = 0.99; // set by model

            soil_model->reset(soil_data);
        }

    }

    void update_water_content(LocalElementAccessorBase<3> ele) override {
        reset_soil_model(ele);
        double storativity = this->ad_->storativity.value(ele.centre(), ele.element_accessor());
        for (unsigned int i=0; i<ele.element_accessor()->n_sides(); i++) {
            double capacity = 0;
            double water_content = 0;
            double phead = ad_->phead_edge_[ele.edge_local_idx(i)];
            if (genuchten_on) {

                  fadbad::B<double> x_phead(phead);
                  fadbad::B<double> evaluated( soil_model->water_content_diff(x_phead) );
                  evaluated.diff(0,1);
                  water_content = evaluated.val();
                  capacity = x_phead.d(0);
            }
            ad_->capacity[ele.side_local_idx(i)] = capacity + storativity;
            ad_->water_content_previous_it[ele.side_local_idx(i)] = water_content + storativity * phead;
        }
    }

    void assemble_sides(LocalElementAccessorBase<3> ele) override
    {
        reset_soil_model(ele);
        cross_section = this->ad_->cross_section.value(ele.centre(), ele.element_accessor());


        double conductivity, head;
        if (genuchten_on) {
            conductivity=0;
            head=0;
            for (unsigned int i=0; i<ele.element_accessor()->n_sides(); i++)
            {
                uint local_edge = ele.edge_local_idx(i);
                double phead = ad_->phead_edge_[local_edge];
                conductivity += soil_model->conductivity(phead);
                head += ad_->phead_edge_[local_edge];
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

        // set lumped source
        double diagonal_coef = ele.measure() * cross_section / ele.n_sides();


        double source_diagonal = diagonal_coef * this->ad_->water_source_density.value(ele.centre(), ele.element_accessor());

        update_water_content(ele);
        for (unsigned int i=0; i<ele.element_accessor()->n_sides(); i++)
        {

            uint local_edge = ele.edge_local_idx(i);
            uint local_side = ele.side_local_idx(i);
            if (this->dirichlet_edge[i] == 0) {

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


                double mass_rhs = mass_diagonal * ad_->phead_edge_[local_edge] / this->ad_->time_step_
                                  + diagonal_coef * water_content_diff / this->ad_->time_step_;

//                 DBGCOUT(<< "source [" << loc_system_.row_dofs[this->loc_edge_dofs[i]] << ", " << loc_system_.row_dofs[this->loc_edge_dofs[i]]
//                             << "] mat: " << -mass_diagonal/this->ad_->time_step_
//                             << " rhs: " << -source_diagonal - mass_rhs
//                             << "\n");
                this->loc_system_.add_value(this->loc_edge_dofs[i], this->loc_edge_dofs[i],
                                            -mass_diagonal/this->ad_->time_step_,
                                            -source_diagonal - mass_rhs);
            }

            if (ad_->balance != nullptr) {
                ad_->balance->add_mass_vec_value(ad_->water_balance_idx, ele.region().bulk_idx(),
                        diagonal_coef*ad_->water_content_previous_it[local_side]);
                ad_->balance->add_source_values(ad_->water_balance_idx, ele.region().bulk_idx(), {(LongIdx)local_edge},{0},{source_diagonal});
            }
        }

    }

    AssemblyDataPtr ad_;

    bool genuchten_on;
    double cross_section;
    std::shared_ptr<SoilModelBase> soil_model;
};



#endif /* SRC_FLOW_ASSEMBLY_LMH_HH_ */
