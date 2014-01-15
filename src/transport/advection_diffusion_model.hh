/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief Discontinuous Galerkin method for equation of transport with dispersion.
 *  @author Jan Stebel
 */

#ifndef AD_MODEL_HH_
#define AD_MODEL_HH_

#include <armadillo>
#include <vector>
#include "input/input_type.hh"


class AdvectionDiffusionModel {
public:

	virtual void init_data(Mesh *mesh_, unsigned int n_subst_, Input::Array bulk_list, Input::Array bc_list) = 0;

//	virtual void set_time(TimeGovernor &time) = 0;

	virtual bool mass_matrix_changed() = 0;

	virtual bool stiffness_matrix_changed() = 0;

	virtual bool rhs_changed() = 0;

	virtual void compute_mass_matrix_coefficient(const std::vector<arma::vec3 > &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<double> &mm_coef) = 0;

	virtual void compute_advection_diffusion_coefficients(const std::vector<arma::vec3> &point_list,
			const std::vector<arma::vec3> &velocity,
			const ElementAccessor<3> &ele_acc,
			std::vector<std::vector<arma::vec3> > &ad_coef,
			std::vector<std::vector<arma::mat33> > &dif_coef) = 0;

	virtual void compute_source_coefficients(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<arma::vec> &sources_conc,
			std::vector<arma::vec> &sources_density,
			std::vector<arma::vec> &sources_sigma) = 0;

	virtual void compute_sources_sigma(const std::vector<arma::vec3> &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<arma::vec> &sources_sigma) = 0;

	virtual ~AdvectionDiffusionModel() {};



};

namespace IT=Input::Type;

class ConcentrationTransportModel : public AdvectionDiffusionModel {
public:

	class ModelEqData : public TransportBase::TransportEqData {
	public:

        enum BC_Type {
            inflow=0,
            dirichlet=1,
            neumann=2,
            robin=3
        };
        static Input::Type::Selection bc_type_selection;


		Field<3, FieldValue<3>::Vector> disp_l;     ///< Longitudal dispersivity (for each substance).
		Field<3, FieldValue<3>::Vector> disp_t;     ///< Transversal dispersivity (for each substance).
		Field<3, FieldValue<3>::Vector> diff_m;     ///< Molecular diffusivity (for each substance).
//		Field<3, FieldValue<3>::Vector> sigma_c;    ///< Transition parameter for diffusive transfer on fractures (for each substance).
//		Field<3, FieldValue<3>::Vector> dg_penalty; ///< Penalty enforcing inter-element continuity of solution (for each substance).

//        BCField<3, FieldValue<3>::EnumVector > bc_type;
//        BCField<3, FieldValue<3>::Vector > bc_flux;
//        BCField<3, FieldValue<3>::Vector > bc_robin_sigma;

        ModelEqData() : TransportBase::TransportEqData("TransportDG")
        {
        	ADD_FIELD(disp_l, "Longitudal dispersivity (for each substance).", IT::Default("0"));
        	ADD_FIELD(disp_t, "Transversal dispersivity (for each substance).", IT::Default("0"));
        	ADD_FIELD(diff_m, "Molecular diffusivity (for each substance).", IT::Default("0"));
//        	ADD_FIELD(sigma_c, "Coefficient of diffusive transfer through fractures (for each substance).", IT::Default("0"));
//        	ADD_FIELD(dg_penalty, "Penalty parameter influencing the discontinuity of the solution (for each substance). "
//        			"Its default value 1 is sufficient in most cases. Higher value diminishes the inter-element jumps.", IT::Default("1.0"));

//            ADD_FIELD(bc_type,"Boundary condition type, possible values: inflow, dirichlet, neumann, robin.", IT::Default("inflow") );
//            bc_type.set_selection(&bc_type_selection);
//            ADD_FIELD(bc_flux,"Flux in Neumann boundary condition.", IT::Default("0.0"));
//            ADD_FIELD(bc_robin_sigma,"Conductivity coefficient in Robin boundary condition.", IT::Default("0.0"));
        }

        RegionSet read_boundary_list_item(Input::Record rec) {
        	// Base method EqDataBase::read_boundary_list_item must be called first!
        	RegionSet domain = EqDataBase::read_boundary_list_item(rec);
            FilePath bcd_file;

            // read transport boundary conditions using old file format .tbc
            if (rec.opt_val("old_boundary_file", bcd_file) )
                OldBcdInput::instance()->read_transport(bcd_file, bc_conc);

            return domain;
        }

        const std::vector<FieldCommonBase *> get_field_list() {
        	return field_list;
        }

	};

protected:

	virtual ModelEqData &data() = 0;

	bool flux_changed;


public:

	ConcentrationTransportModel() {};

	void init_data(Mesh *mesh_, unsigned int n_subst_, Input::Array bulk_list, Input::Array bc_list) {
//		data().set_mesh(mesh_);
		data().init_conc.set_n_comp(n_subst_);
		data().bc_conc.set_n_comp(n_subst_);
	    data().sources_density.set_n_comp(n_subst_);
	    data().sources_sigma.set_n_comp(n_subst_);
	    data().sources_conc.set_n_comp(n_subst_);
		data().diff_m.set_n_comp(n_subst_);
		data().disp_l.set_n_comp(n_subst_);
		data().disp_t.set_n_comp(n_subst_);
//		data().init_from_input(bulk_list, bc_list);
		flux_changed = true;
	}

//	void set_time(TimeGovernor &time) { data().set_time(time); }

	bool mass_matrix_changed() {
		return (data().cross_section->changed() || data().por_m.changed());
	}

	bool stiffness_matrix_changed() {
		return (flux_changed ||
				data().disp_l.changed() ||
				data().disp_t.changed() ||
				data().diff_m.changed() ||
				data().por_m.changed() ||
				data().cross_section->changed());
	}

	bool rhs_changed() {
		return (flux_changed ||
				data().bc_conc.changed() ||
				data().sources_conc.changed() ||
				data().sources_density.changed() ||
				data().sources_sigma.changed());
	}

	void compute_mass_matrix_coefficient(const std::vector<arma::vec3 > &point_list,
			const ElementAccessor<3> &ele_acc,
			std::vector<double> &mm_coef)
	{
		vector<double> elem_csec(point_list.size()), por_m(point_list.size());

		data().cross_section->value_list(point_list, ele_acc, elem_csec);
		data().por_m.value_list(point_list, ele_acc, por_m);

		for (unsigned int i=0; i<mm_coef.size(); i++)
			mm_coef[i] = elem_csec[i]*por_m[i];
	}

	void calculate_dispersivity_tensor(arma::mat33 &K, const arma::vec3 &velocity,
			double Dm, double alphaL, double alphaT, double porosity, double cross_cut)
	{
	    double vnorm = arma::norm(velocity, 2);

		if (fabs(vnorm) > sqrt(numeric_limits<double>::epsilon()))
			for (int i=0; i<3; i++)
				for (int j=0; j<3; j++)
					K(i,j) = velocity[i]*velocity[j]/(vnorm*vnorm)*(alphaL-alphaT) + alphaT*(i==j?1:0);
		else
			K.zeros();

		K = ((vnorm*porosity*cross_cut)*K + (Dm*pow(porosity, 1./3)*porosity*cross_cut)*arma::eye(3,3));
	}

	void compute_advection_diffusion_coefficients(const std::vector<arma::vec3 > &point_list,
			const std::vector<arma::vec3> &velocity,
			const ElementAccessor<3> &ele_acc,
			std::vector<std::vector<arma::vec3> > &ad_coef,
			std::vector<std::vector<arma::mat33> > &dif_coef)
	{
		const unsigned int qsize = point_list.size();
		const unsigned int n_subst = dif_coef.size();
		std::vector<arma::vec> Dm(qsize), alphaL(qsize), alphaT(qsize);
		std::vector<double> por_m(qsize), csection(qsize);

		data().diff_m.value_list(point_list, ele_acc, Dm);
		data().disp_l.value_list(point_list, ele_acc, alphaL);
		data().disp_t.value_list(point_list, ele_acc, alphaT);
		data().por_m.value_list(point_list, ele_acc, por_m);
		data().cross_section->value_list(point_list, ele_acc, csection);

		for (unsigned int i=0; i<qsize; i++) {
			for (int sbi=0; sbi<n_subst; sbi++) {
				ad_coef[sbi][i] = velocity[i];
				calculate_dispersivity_tensor(dif_coef[sbi][i], velocity[i], Dm[i][sbi], alphaL[i][sbi], alphaT[i][sbi], por_m[i], csection[i]);
			}
		}
	}

	void compute_source_coefficients(const std::vector<arma::vec3> &point_list,
				const ElementAccessor<3> &ele_acc,
				std::vector<arma::vec> &sources_conc,
				std::vector<arma::vec> &sources_density,
				std::vector<arma::vec> &sources_sigma)
	{
		data().sources_conc.value_list(point_list, ele_acc, sources_conc);
		data().sources_density.value_list(point_list, ele_acc, sources_density);
		data().sources_sigma.value_list(point_list, ele_acc, sources_sigma);
	}

	void compute_sources_sigma(const std::vector<arma::vec3> &point_list,
				const ElementAccessor<3> &ele_acc,
				std::vector<arma::vec> &sources_sigma)
	{
		data().sources_sigma.value_list(point_list, ele_acc, sources_sigma);
	}

	~ConcentrationTransportModel() {};
};



#endif /* AD_MODEL_HH_ */
