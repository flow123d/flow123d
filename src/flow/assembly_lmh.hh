/*
 * assembly_lmh.hh
 *
 *  Created on: Apr 30, 2016
 *      Author: jb
 */

#ifndef SRC_FLOW_ASSEMBLY_LMH_HH_
#define SRC_FLOW_ASSEMBLY_LMH_HH_

#include "flow/darcy_flow_assembly.hh"
#include "soil_models.hh"




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
      system_(data->system_)
    {}

    void reset_soil_model(LocalElementAccessor ele) {
        genuchten_on = (this->ad_->genuchten_p_head_scale.field_result({ele.iter->region()}) != result_zeros);
        //DBGMSG("g on: %d\n", genuchten_on);
        if (genuchten_on) {
            SoilData soil_data = {
                    this->ad_->genuchten_n_exponent.value(ele.centre(), ele.accessor()),
                    this->ad_->genuchten_p_head_scale.value(ele.centre(), ele.accessor()),
                    this->ad_->water_content_residual.value(ele.centre(), ele.accessor()),
                    this->ad_->water_content_saturated.value(ele.centre(), ele.accessor()),
                    this->ad_->conductivity.value(ele.centre(), ele.accessor()),
                    1.001, // cut fraction, todo: fix meaning
            };

            soil_model.reset(soil_data);
        }

    }

    void assembly_local_matrix(LocalElementAccessor ele) override
    {
        ele.update();
        reset_soil_model(ele);
        cross_section = this->ad_->cross_section.value(ele.centre(), ele.accessor());


        double conductivity;
        if (genuchten_on) {
            conductivity=0;
            FOR_ELEMENT_SIDES(ele.iter, i)
            {
                uint local_edge = ele.local_edge_idx[i];
                conductivity += soil_model.conductivity(ad_->phead_edge_[local_edge]);
            }
        } else {
            conductivity = this->ad_->conductivity.value(ele.centre(), ele.accessor());
        }

        double scale = 1 / cross_section / conductivity;
        *(system_.local_matrix) = scale * this->assembly_local_geometry_matrix(ele.iter);

        assembly_source_term(ele);
    }

    /***
     * Called from assembly_local_matrix, assumes precomputed:
     * cross_section, genuchten_on, soil_model
     */
    void assembly_source_term(LocalElementAccessor ele) {

        // set lumped source

        double storativity = this->ad_->storativity.value(ele.centre(), ele.accessor());
        double diagonal_coef = ele.iter->measure() * cross_section / ele.iter->n_sides();
        double capacity = 0;
        double water_content_diff = 0;
        double source_diagonal = diagonal_coef * this->ad_->water_source_density.value(ele.centre(), ele.accessor());

        FOR_ELEMENT_SIDES(ele.iter, i)
        {
            uint local_edge = ele.local_edge_idx[i];
            uint edge_row = ele.edge_row[i];

            if (genuchten_on) {
                fadbad::B<double> x_phead(ad_->phead_edge_[local_edge]);
                fadbad::B<double> evaluated( soil_model.water_content(x_phead) );
                evaluated.diff(0,1);
                capacity =  x_phead.d(0);

                ad_->water_content_previous_it[local_edge] = diagonal_coef * evaluated.val();
                water_content_diff = -ad_->water_content_previous_it[local_edge] + ad_->water_content_previous_time[local_edge];
            }

            double mass_balance_diagonal = diagonal_coef * (capacity + storativity);
            double mass_diagonal = mass_balance_diagonal / this->ad_->time_step_;

            //cout << "mesh edge: " << mesh_edge << "local: " << local_edge << endl;
            double mass_rhs = mass_diagonal * ad_->phead_edge_[local_edge] + water_content_diff / this->ad_->time_step_;

            system_.lin_sys->mat_set_value(edge_row, edge_row, -mass_diagonal );
            system_.lin_sys->rhs_set_value(edge_row, -source_diagonal - mass_rhs);

            if (system_.balance != nullptr) {
                system_.balance->add_mass_matrix_values(ad_->water_balance_idx_, ele.iter->region().bulk_idx(), {edge_row}, {mass_balance_diagonal});
                system_.balance->add_source_rhs_values(ad_->water_balance_idx_, ele.iter->region().bulk_idx(), {edge_row}, {source_diagonal});
            }
        }

    }

    AssemblyDataPtr ad_;
    RichardsSystem system_;

    SoilModel_VanGenuchten soil_model;
    bool genuchten_on;
    double cross_section;
    /***
     * Called from assembly_local_matrix, assumes precomputed:
     * cross_section, genuchten_on, soil_model
     */
    void assembly_source_term(LocalElementAccessor ele) {

        // set lumped source

        double storativity = this->ad_->storativity.value(ele.centre(), ele.accessor());
        double diagonal_coef = ele.iter->measure() * cross_section / ele.iter->n_sides();
        double capacity = 0;
        double water_content_diff = 0;
        double source_diagonal = diagonal_coef * this->ad_->water_source_density.value(ele.centre(), ele.accessor());
        double side_water_content=0;

        FOR_ELEMENT_SIDES(ele.iter, i)
        {
            uint local_edge = ele.local_edge_idx[i];
            uint local_side = ele.local_side_idx[i];
            uint edge_row = ele.edge_row[i];
            double water_content=0;

            if (genuchten_on) {
                fadbad::B<double> x_phead(ad_->phead_edge_[local_edge]);
                fadbad::B<double> evaluated( soil_model.water_content(x_phead) );
                evaluated.diff(0,1);
                capacity =  x_phead.d(0);

                water_content = evaluated.val();
            }
            side_water_content = diagonal_coef * (water_content + storativity*ad_->phead_edge_[local_edge]);


            //ad_->water_content_previous_it[local_side] = side_water_content;


            water_content_diff = -side_water_content + ad_->water_content_previous_time[local_side];

            double mass_balance_diagonal = diagonal_coef * (capacity + storativity);
            double mass_diagonal = mass_balance_diagonal / this->ad_->time_step_;

            //cout << "mesh edge: " << mesh_edge << "local: " << local_edge << endl;
            double mass_rhs = mass_diagonal * ad_->phead_edge_[local_edge] + water_content_diff / this->ad_->time_step_;
            //DBGMSG("rhs: %d %f\n", local_side, ad_->phead_edge_[local_edge]);


            system_.lin_sys->mat_set_value(edge_row, edge_row, -mass_diagonal );
            system_.lin_sys->rhs_set_value(edge_row, source_diagonal - mass_rhs);

            if (system_.balance != nullptr) {
                system_.balance->add_mass_matrix_values(ad_->water_balance_idx_, ele.iter->region().bulk_idx(), {edge_row}, {mass_balance_diagonal});
                system_.balance->add_source_rhs_values(ad_->water_balance_idx_, ele.iter->region().bulk_idx(), {edge_row}, {source_diagonal});
            }
        }

    }

    void init_water_content(LocalElementAccessor ele, double p_head) {
        ele.update();
        reset_soil_model(ele);
        double storativity = this->ad_->storativity.value(ele.centre(), ele.accessor());
        cross_section = this->ad_->cross_section.value(ele.centre(), ele.accessor());
        double diagonal_coef = ele.iter->measure() * cross_section / ele.iter->n_sides();

        double water_content = 0.0;
        if (genuchten_on) {
            water_content = soil_model.water_content(p_head);
        }

        for(int i=0; i<ele.iter->n_sides(); i++) {
            uint local_side = ele.local_side_idx[i];
            double wc_init = diagonal_coef * (water_content + storativity * p_head);
            ad_->water_content_previous_it[local_side] = wc_init;
            ad_->water_content_previous_time[local_side] = wc_init;
            //DBGMSG("wc_init: %d %f\n", local_side, ad_->water_content_previous_it[local_side]);
        }
    }

    AssemblyDataPtr ad_;
    RichardsSystem system_;

    SoilModel_VanGenuchten soil_model;
    bool genuchten_on;
    double cross_section;
};



#endif /* SRC_FLOW_ASSEMBLY_LMH_HH_ */
