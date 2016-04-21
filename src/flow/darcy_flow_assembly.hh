/*
 * darcy_flow_assembly.hh
 *
 *  Created on: Apr 21, 2016
 *      Author: jb
 */

#ifndef SRC_FLOW_DARCY_FLOW_ASSEMBLY_HH_
#define SRC_FLOW_DARCY_FLOW_ASSEMBLY_HH_

#include <memory>
#include "fem/mapping_p1.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/fe_rt.hh"
#include "quadrature/quadrature_lib.hh"


    struct AssemblyData
    {
        Mesh *mesh;
        DarcyFlowMH_Steady::EqData* data;
        MH_DofHandler *mh_dh;

        double effective_conductivity;
        double cross_section;
        double edge_saturation[4];

        LinSys *system;

    };


    class AssemblyBase
    {
    public:
        typedef std::vector<std::shared_ptr<AssemblyBase> > MultidimAssembly;

        virtual ~AssemblyBase() {}

        /**
         * Generic creator of multidimensional assembly, i.e. vector of
         * particular assembly objects.
         */
        template< template<int dim> class Impl >
        static MultidimAssembly create(Mesh &mesh, DarcyFlowMH_Steady::EqData& data, MH_DofHandler &mh_dh) {
            AssemblyData ad;
            ad.mesh = &mesh;
            ad.data = &data;
            ad.mh_dh = &mh_dh;
            return { std::make_shared<Impl<1> >(ad),
                std::make_shared<Impl<2> >(ad),
                std::make_shared<Impl<3> >(ad) };

        }


        // assembly just A block of local matrix
        virtual void assembly_local_matrix(ElementFullIter ele,
                                           arma::mat &local_matrix) =0;

        // assembly compatible neighbourings
        virtual void assembly_local_vb(double *local_vb,
                                       ElementFullIter ele,
                                       Neighbour *ngh) = 0;

        // compute velocity value in the barycenter
        // TODO: implement and use general interpolations between discrete spaces
        virtual arma::vec3 make_element_vector(ElementFullIter ele) = 0;
    };

    template<int dim>
    class AssemblyMH : public AssemblyBase
    {
    public:
        AssemblyMH<dim>(AssemblyData ad)
        : quad_(3),
          fe_values_(map_, quad_, fe_rt_,
                    update_values | update_gradients | update_JxW_values | update_quadrature_points),

          side_quad_(1),
          fe_side_values_(map_, side_quad_, fe_p_disc_, update_normal_vectors),

          velocity_interpolation_quad_(0), // veloctiy values in barycenter
          velocity_interpolation_fv_(map_,velocity_interpolation_quad_, fe_rt_, update_values | update_quadrature_points),

          ad_(ad)
        {}


        ~AssemblyMH<dim>() override
        {}

        void assembly_local_matrix(ElementFullIter ele,
                arma::mat &local_matrix) override
        {
            //START_TIMER("Assembly<dim>::assembly_local_matrix");
            fe_values_.reinit(ele);
            unsigned int ndofs = fe_values_.get_fe()->n_dofs();
            unsigned int qsize = fe_values_.get_quadrature()->size();
            local_matrix.zeros(ndofs, ndofs);

            double scale = 1
                           / ad_.data->conductivity.value( ele->centre(), ele->element_accessor() )
                           / ad_.data->cross_section.value( ele->centre(), ele->element_accessor() );

            for (unsigned int k=0; k<qsize; k++)
            {
                for (unsigned int i=0; i<ndofs; i++)
                {
                     for (unsigned int j=0; j<ndofs; j++)
                        local_matrix[i*ndofs+j] +=
                                scale
                                * arma::dot(fe_values_.shape_vector(i,k),
                                            (ad_.data->anisotropy.value(ele->centre(), ele->element_accessor() )).i()
                                             * fe_values_.shape_vector(j,k)
                                           )
                                * fe_values_.JxW(k);
                }
            }
        }


        void assembly_local_vb(double *local_vb,  ElementFullIter ele, Neighbour *ngh) override
        {
            //START_TIMER("Assembly<dim>::assembly_local_vb");
            // compute normal vector to side
            arma::vec3 nv;
            ElementFullIter ele_higher = ad_.mesh->element.full_iter(ngh->side()->element());
            fe_side_values_.reinit(ele_higher, ngh->side()->el_idx());
            nv = fe_side_values_.normal_vector(0);

            double value = ad_.data->sigma.value( ele->centre(), ele->element_accessor()) *
                            2*ad_.data->conductivity.value( ele->centre(), ele->element_accessor()) *
                            arma::dot(ad_.data->anisotropy.value( ele->centre(), ele->element_accessor())*nv, nv) *
                            ad_.data->cross_section.value( ngh->side()->centre(), ele_higher->element_accessor() ) * // cross-section of higher dim. (2d)
                            ad_.data->cross_section.value( ngh->side()->centre(), ele_higher->element_accessor() ) /
                            ad_.data->cross_section.value( ele->centre(), ele->element_accessor() ) *      // crossection of lower dim.
                            ngh->side()->measure();

            local_vb[0] = -value;   local_vb[1] = value;
            local_vb[2] = value;    local_vb[3] = -value;
        }


        arma::vec3 make_element_vector(ElementFullIter ele) override
        {
            //START_TIMER("Assembly<dim>::make_element_vector");
            arma::vec3 flux_in_center;
            flux_in_center.zeros();

            velocity_interpolation_fv_.reinit(ele);
            for (unsigned int li = 0; li < ele->n_sides(); li++) {
                flux_in_center += ad_.mh_dh->side_flux( *(ele->side( li ) ) )
                          * velocity_interpolation_fv_.shape_vector(li,0);
            }

            flux_in_center /= ad_.data->cross_section.value(ele->centre(), ele->element_accessor() );
            return flux_in_center;
        }


        // assembly volume integrals
        FE_RT0<dim,3> fe_rt_;
        MappingP1<dim,3> map_;
        QGauss<dim> quad_;
        FEValues<dim,3> fe_values_;

        // assembly face integrals (BC)
        QGauss<dim-1> side_quad_;
        FE_P_disc<0,dim,3> fe_p_disc_;
        FESideValues<dim,3> fe_side_values_;

        // Interpolation of velocity into barycenters
        QGauss<dim> velocity_interpolation_quad_;
        FEValues<dim,3> velocity_interpolation_fv_;

        // data shared by assemblers of different dimension
        AssemblyData ad_;

    };



    template<int dim>
    class AssemblyLMH : public AssemblyMH<dim> {

    };


#endif /* SRC_FLOW_DARCY_FLOW_ASSEMBLY_HH_ */
