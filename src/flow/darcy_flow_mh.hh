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
 * @brief mixed-hybrid model of linear Darcy flow, possibly unsteady.
 *
 * Main object for mixed-hybrid discretization of the linear elliptic PDE (Laplace)
 * on  a multidimensional domain. Discretization of saturated Darcy flow using
 * RT0 approximation for the velocity
 *
 *
 *
 *  @author Jan Brezina
 *
 */

/*
 * list of files dependent on this one:
 *
 * posprocess.cc
 * problem.cc
 * main.hh
 * transport.cc
 */


#ifndef DARCY_FLOW_MH_HH
#define DARCY_FLOW_MH_HH

#include <memory>
#include "input/input_type.hh"

#include <petscmat.h>
#include "system/sys_vector.hh"
#include "coupling/equation.hh"
#include "flow/mh_dofhandler.hh"
#include "input/input_type.hh"
#include "la/linsys_BDDC.hh"
#include "la/linsys_PETSC.hh"

#include "fields/bc_field.hh"
#include "fields/field.hh"
#include "fields/field_set.hh"
#include "fields/field_add_potential.hh"
#include "flow/old_bcd.hh"

class A;

/// external types:
class LinSys;
struct Solver;
class Mesh;
class SchurComplement;
class Distribution;
class SparseGraph;
class LocalToGlobalMap;
class DarcyFlowMHOutput;
class Balance;
class VectorSeqDouble;

template<unsigned int dim, unsigned int spacedim> class FE_RT0;
template<unsigned int degree, unsigned int dim, unsigned int spacedim> class FE_P_disc;
template<unsigned int dim, unsigned int spacedim> class MappingP1;
template<unsigned int dim, unsigned int spacedim> class FEValues;

/**
 * @brief Mixed-hybrid model of linear Darcy flow, possibly unsteady.
 *
 * Abstract class for various implementations of Darcy flow. In future there should be
 * one further level of abstraction for general time dependent problem.
 *
 * maybe TODO:
 * split compute_one_step to :
 * 1) prepare_next_timestep
 * 2) actualize_solution - this is for iterative nonlinear solvers
 *
 */

class DarcyFlowMH : public EquationBase {
public:
    enum MortarMethod {
        NoMortar = 0,
        MortarP0 = 1,
        MortarP1 = 2
    };
    
    /** @brief Data for Darcy flow equation.
     *  
     */
    class EqData : public FieldSet {
    public:

        /**
         * For compatibility with old BCD file we have to assign integer codes starting from 1.
         */
        enum BC_Type {
            none=0,
            dirichlet=1,
            neumann=2,
            robin=3,
            total_flux=4
        };
        static const Input::Type::Selection & get_bc_type_selection();

        /// Collect all fields
        EqData();


        Field<3, FieldValue<3>::TensorFixed > anisotropy;
        Field<3, FieldValue<3>::Scalar > conductivity;
        Field<3, FieldValue<3>::Scalar > cross_section;
        Field<3, FieldValue<3>::Scalar > water_source_density;
        Field<3, FieldValue<3>::Scalar > sigma;

        BCField<3, FieldValue<3>::Enum > bc_type; // Discrete need Selection for initialization
        BCField<3, FieldValue<3>::Scalar > bc_pressure; 
        BCField<3, FieldValue<3>::Scalar > bc_flux;
        BCField<3, FieldValue<3>::Scalar > bc_robin_sigma;
        
        //TODO: these belong to Unsteady flow classes
        //as long as Unsteady is descendant from Steady, these cannot be transfered..
        Field<3, FieldValue<3>::Scalar > init_pressure;
        Field<3, FieldValue<3>::Scalar > storativity;

        /**
         * Gravity vector and constant shift of pressure potential. Used to convert piezometric head
         * to pressure head and vice versa.
         */
        arma::vec4 gravity_;

        FieldSet	time_term_fields;
        FieldSet	main_matrix_fields;
        FieldSet	rhs_fields;
    };


    /**
     * Model for transition coefficients due to Martin, Jaffre, Roberts (see manual for full reference)
     *
     * TODO:
     * - how we can reuse values computed during assembly
     *   we want to make this class see values in
     *
     */
    DarcyFlowMH(Mesh &mesh, const Input::Record in_rec)
    : EquationBase(mesh, in_rec)
    {}

    static const Input::Type::Selection & get_mh_mortar_selection();
    static Input::Type::AbstractRecord & get_input_type();

    void get_velocity_seq_vector(Vec &velocity_vec)
        { velocity_vec = velocity_vector; }

    const MH_DofHandler &get_mh_dofhandler() {
        double *array;
        unsigned int size;
        get_solution_vector(array, size);

        // here assume that velocity field is extended as constant
        // to the previous time, so here we set left bound of the interval where the velocity
        // has current value; this may not be good for every transport !!
        // we can resolve this when we use FieldFE to store computed velocities in few last steps and
        // let every equation set time according to nature of the time scheme

        // in particular this setting is necessary to prevent ConvectinTransport to recreate the transport matrix
        // every timestep ( this may happen for unsteady flow if we would use time->t() here since it returns infinity.
        mh_dh.set_solution(time_->last_t(), array, solution_precision());
       return mh_dh;
    }
    
    virtual void set_concentration_vector(Vec &vc){};


protected:
    void setup_velocity_vector() {
        double *velocity_array;
        unsigned int size;

        get_solution_vector(velocity_array, size);
        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, mesh_->n_sides(), velocity_array, &velocity_vector);

    }

    virtual double solution_precision() const = 0;

    bool solution_changed_for_scatter;
    Vec velocity_vector;
    MH_DofHandler mh_dh;    // provides access to seq. solution fluxes and pressures on sides

    MortarMethod mortar_method_;

    /// object for calculation and writing the water balance to file.
    boost::shared_ptr<Balance> balance_;
    /// index of water balance within the Balance object.
    unsigned int water_balance_idx_;
};


/**
 * @brief Mixed-hybrid of steady Darcy flow with sources and variable density.
 *
 * solve equations:
 * @f[
 *      q= -{\mathbf{K}} \nabla h -{\mathbf{K}} R \nabla z
 * @f]
 * @f[
 *      \mathrm{div} q = f
 * @f]
 *
 * where
 * - @f$ q @f$ is flux @f$[ms^{-1}]@f$ for 3d, @f$[m^2s^{-1}]@f$ for 2d and @f$[m^3s^{-1}]@f$ for 1d.
 * - @f$ \mathbf{K} @f$ is hydraulic tensor ( its orientation for 2d, 1d case is questionable )
 * - @f$ h = \frac{\pi}{\rho_0 g}+z @f$ is pressure head, @f$ \pi, \rho_0, g @f$ are the pressure, water density, and acceleration of gravity , respectively.
 *   Assumes gravity force acts counter to the direction of the @f$ z @f$ axis.
 * - @f$ R @f$ is destity or gravity variability coefficient. For density driven flow it should be
 * @f[
 *    R = \frac{\rho}{\rho_0} -1 = \rho_0^{-1}\sum_{i=1}^s c_i
 * @f]
 *   where @f$ c_i @f$ is concentration in @f$ kg m^{-3} @f$.
 *
 *
 *   TODO:
 *   - consider create auxiliary classes for creation of inverse of block A
 */
class DarcyFlowMH_Steady : public DarcyFlowMH
{
public:
	typedef DarcyFlowMH FactoryBaseType;
  
    class EqData : public DarcyFlowMH::EqData {
    public:
      
      EqData() : DarcyFlowMH::EqData()
      {}
    };
    
    DarcyFlowMH_Steady(Mesh &mesh, const Input::Record in_rec, bool make_tg=true);

    static const Input::Type::Record & get_input_type();

    virtual void update_solution();
    virtual void get_solution_vector(double * &vec, unsigned int &vec_size);
    virtual void get_parallel_solution_vector(Vec &vector);
    
    /// postprocess velocity field (add sources)
    virtual void postprocess();
    virtual void output_data() override;

    ~DarcyFlowMH_Steady();


protected:
    class AssemblyBase;
    template<unsigned int dim> class Assembly;
    
    class AssemblyData
    {
    public:
        AssemblyData(Mesh *mesh,
                     EqData *data,
                     LinSys *ls, 
                     Distribution *edge_dist, 
                     Distribution *el_dist,
                     Distribution *side_dist,
                     boost::shared_ptr<Balance> balance,
                     MH_DofHandler *mh_dh,
                     unsigned int water_balance_idx,
                     int n_schur,
                     int *el_for_loc,
                     int *row_4_el,
                     int *side_row_4_id,
                     int *row_4_edge
                     );
    private:
        Mesh *mesh;
        LinSys *ls;
        Distribution *edge_ds;          
        Distribution *el_ds;            
        Distribution *side_ds;          
        EqData* data;
        boost::shared_ptr<Balance> balance;
        unsigned int water_balance_idx;
        MH_DofHandler *mh_dh;
        
        int n_schur_compls;
        int *el_4_loc;
        int *row_4_el;
        int *side_row_4_id;
        int *row_4_edge;
        
    template<unsigned int dim>
    friend class Assembly;
    };
    
    class AssemblyBase
    {
    public:
        virtual void assembly_matrix() = 0;
        virtual void make_element_vector(VectorSeqDouble &ele_flux) = 0;
        void set_data(AssemblyData *data);
    protected:
        AssemblyData *d;
    };
    
    template<unsigned int dim>
    class Assembly : public AssemblyBase
    {
    public:
        Assembly<dim>();
        ~Assembly<dim>();
        void assembly_matrix() override;
        void make_element_vector(VectorSeqDouble &ele_flux) override;
        void assembly_local_matrix(arma::mat &local_matrix, 
                                   ElementFullIter ele, 
                                   FEValues<dim,3> & fe_values);
    private:
        FE_RT0<dim,3> *fe_rt_;
        MappingP1<dim,3> *map_;
    };
    
    void make_serial_scatter();
    virtual void modify_system()
    { ASSERT(0, "Modify system called for Steady darcy.\n"); };
    virtual void setup_time_term()
    { ASSERT(0, "Setup time term called for Steady darcy.\n"); };


    void prepare_parallel( const Input::AbstractRecord in_rec);
    void make_row_numberings();

    /**
     * Create and preallocate MH linear system (including matrix, rhs and solution vectors)
     */
    void create_linear_system();

    /**
     * Read initial condition into solution vector.
     * Must be called after create_linear_system.
     *
     */
    virtual void read_init_condition() {};

    /**
     * Abstract assembly method used for both assembly and preallocation.
     * Assembly only steady part of the equation.
     * TODO:
     * - use general preallocation methods in DofHandler
     * - include alos time term
     * - add support for Robin type sources
     * - support for nonlinear solvers - assembly either residual vector, matrix, or both (using FADBAD++)
     *
     */
    void assembly_steady_mh_matrix();
    
    /** Assembly of a local mass matrix on and element.
     * Auxiliary function for @p assembly_steady_mh_matrix.
     */
    template<unsigned int dim>
    void assembly_steady_mh_local_matrix(arma::mat &local_matrix, ElementFullIter ele, 
                                         FEValues<dim,3> & fe_values);

    /// Source term is implemented differently in LMH version.
    virtual void assembly_source_term();

    /**
     * Assembly or update whole linear system.
     */
    void assembly_linear_system();

    void set_mesh_data_for_bddc(LinSys_BDDC * bddc_ls);
    double solution_precision() const;


    DarcyFlowMHOutput *output_object;

	int size;				// global size of MH matrix
	int  n_schur_compls;  	// number of shur complements to make
	double  *solution; 			// sequantial scattered solution vector


	LinSys *schur0;  		//< whole MH Linear System

	AssemblyData *assembly_data_;
	std::vector<AssemblyBase *> assembly_;
	
	// parallel
	Distribution *edge_ds;          //< optimal distribution of edges
	Distribution *el_ds;            //< optimal distribution of elements
	Distribution *side_ds;          //< optimal distribution of elements
	boost::shared_ptr<Distribution> rows_ds;          //< final distribution of rows of MH matrix

	int *el_4_loc;		        //< array of idexes of local elements (in ordering matching the optimal global)
	int *row_4_el;		        //< element index to matrix row
	int *side_id_4_loc;		//< array of ids of local sides
	int	*side_row_4_id;		//< side id to matrix row
	int *edge_4_loc;		//< array of indexes of local edges
	int	*row_4_edge;		//< edge index to matrix row

	// MATIS related arrays
        boost::shared_ptr<LocalToGlobalMap> global_row_4_sub_row;           //< global dof index for subdomain index

	// gather of the solution
	Vec sol_vec;			                 //< vector over solution array
	VecScatter par_to_all;
        
  EqData data_;

  friend class DarcyFlowMHOutput;
  friend class P0_CouplingAssembler;
  friend class P1_CouplingAssembler;

private:
  /// Registrar of class to factory
  static const int registrar;
};


class P0_CouplingAssembler {
public:
	P0_CouplingAssembler(const DarcyFlowMH_Steady &darcy)
	: darcy_(darcy),
	  master_list_(darcy.mesh_->master_elements),
	  intersections_(darcy.mesh_->intersections),
	  master_(nullptr),
	  tensor_average(2)
	{
		arma::mat master_map(1,2), slave_map(1,3);
		master_map.fill(1.0 / 2);
		slave_map.fill(1.0 / 3);

		tensor_average[0].push_back( trans( master_map ) * master_map );
		tensor_average[0].push_back( trans( master_map ) * slave_map );
		tensor_average[1].push_back( trans( slave_map ) * master_map );
		tensor_average[1].push_back( trans( slave_map ) * slave_map );
	}

	void assembly(LinSys &ls);
	void pressure_diff(int i_ele,
			vector<int> &dofs,
			unsigned int &ele_type,
			double &delta,
			arma::vec &dirichlet);
private:
	typedef vector<unsigned int> IsecList;

	const DarcyFlowMH_Steady &darcy_;

	const vector<IsecList> &master_list_;
	const vector<Intersection> &intersections_;

	vector<IsecList>::const_iterator ml_it_;
	const Element *master_;

	/// Row matrices to compute element pressure as average of boundary pressures
	vector< vector< arma::mat > > tensor_average;
	/// measure of master element, should be sum of intersection measures
	double delta_0;
};



class P1_CouplingAssembler {
public:
	P1_CouplingAssembler(const DarcyFlowMH_Steady &darcy)
	: darcy_(darcy),
	  intersections_(darcy.mesh_->intersections),
	  rhs(5),
	  dofs(5),
	  dirichlet(5)
	{
		rhs.zeros();
	}

	void assembly(LinSys &ls);
	void add_sides(const Element * ele, unsigned int shift, vector<int> &dofs, vector<double> &dirichlet);
private:

	const DarcyFlowMH_Steady &darcy_;
	const vector<Intersection> &intersections_;

	arma::vec rhs;
	vector<int> dofs;
	vector<double> dirichlet;
};



void mat_count_off_proc_values(Mat m, Vec v);



/**
 * @brief Mixed-hybrid solution of unsteady Darcy flow.
 *
 * Standard discretization with time term and sources picewise constant
 * on the element. This leads to violation of the discrete maximum principle for
 * non-acute meshes or to too small timesteps. For simplicial meshes this can be solved by lumping to the edges. See DarcyFlowLMH_Unsteady.
 */

class DarcyFlowMH_Unsteady : public DarcyFlowMH_Steady
{
public:
	typedef DarcyFlowMH FactoryBaseType;
  
    DarcyFlowMH_Unsteady(Mesh &mesh, const Input::Record in_rec);
    DarcyFlowMH_Unsteady();

    static const Input::Type::Record & get_input_type();
protected:
    void read_init_condition() override;
    void modify_system() override;
    void setup_time_term();
    
private:
    /// Registrar of class to factory
    static const int registrar;

    Vec steady_diagonal;
    Vec steady_rhs;
    Vec new_diagonal;
    Vec previous_solution;

};

/**
 * @brief Edge lumped mixed-hybrid solution of unsteady Darcy flow.
 *
 * The time term and sources are evenly distributed form an element to its edges.
 * This applies directly to the second Schur complement. After this system for pressure traces is solved we reconstruct pressures and side flows as follows:
 *
 * -# Element pressure is  average of edge pressure. This is in fact same as the MH for steady case so we let SchurComplement class do its job.
 *
 * -# We let SchurComplement to reconstruct fluxes and then account time term and sources which are evenly distributed from an element to its sides.
 *    It can be proved, that this keeps continuity of the fluxes over the edges.
 *
 * This lumping technique preserves discrete maximum principle for any time step provided one use acute mesh. But in practice even worse meshes are tractable.
 */
class DarcyFlowLMH_Unsteady : public DarcyFlowMH_Steady
{
public:
	typedef DarcyFlowMH FactoryBaseType;
  
    DarcyFlowLMH_Unsteady(Mesh &mesh, const Input::Record in_rec);
    DarcyFlowLMH_Unsteady();
    
    static const Input::Type::Record & get_input_type();
protected:
    void read_init_condition() override;
    void modify_system() override;
    void assembly_source_term() override;
    void setup_time_term();
    virtual void postprocess();
private:
    /// Registrar of class to factory
    static const int registrar;

    Vec steady_diagonal;
    Vec steady_rhs;
    Vec new_diagonal;
    Vec previous_solution;
    //Vec time_term;
};

#endif  //DARCY_FLOW_MH_HH
//-----------------------------------------------------------------------------
// vim: set cindent:

