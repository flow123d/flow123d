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

    template<int dim>
    class AssemblyLMH : public AssemblyMH<dim> {
    public:
        typedef std::shared_ptr<DarcyFlowLMH_Unsteady::EqData> AssemblyDataPtr;
        AssemblyLMH(AssemblyDataPtr data)
        : AssemblyMH<dim>(data)
        {}

        void assembly_local_matrix(ElementFullIter ele, unsigned int loc_ele_idx,
                                                   arma::mat &local_matrix) override
        {
            double cs = this->ad_->cross_section.value(ele->centre(), ele->element_accessor());

            /*
            SoilData soil_data = {
                    this->ad_.data->genuchten_n_exponent.value(ele->centre(), ele->element_accessor()),
                    this->ad_.data->genuchten_p_head_scale.value(ele->centre(), ele->element_accessor()),
                    this->ad_.data->water_content_residual.value(ele->centre(), ele->element_accessor()),
                    this->ad_.data->water_content_saturated.value(ele->centre(), ele->element_accessor()),
                    1.001, // cut fraction, tod: fix meaning
                    this->ad_.data->conductivity.value(ele->centre(), ele->element_accessor())
            };
            SoilModel_VanGenuchten soil_model(soil_data);
            double conductivity_mean=0;
            FOR_ELEMENT_SIDES(ele,i)
            {
                int mesh_edge=ele->side(i)->edge_idx();
                int local_edge = ad_.data->edge_new_local_4_mesh_idx_[mesh_edge];

                conductivity_mean += soil_model.conductivity(ad_.data->phead_edge_[local_edge]);
            }*/

            double conductivity_mean = this->ad_->conductivity.value(ele->centre(), ele->element_accessor());
            double scale = 1 / cs / conductivity_mean;
            local_matrix = scale * this->assembly_local_geometry_matrix(ele);
        }

    };



#endif /* SRC_FLOW_ASSEMBLY_LMH_HH_ */
