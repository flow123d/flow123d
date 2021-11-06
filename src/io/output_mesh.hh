/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    output_mesh.hh
 * @brief   Classes for auxiliary output mesh.
 */

#ifndef OUTPUT_MESH_HH_
#define OUTPUT_MESH_HH_

#include <memory>                     // for shared_ptr, enable_shared_from_...
#include <string>                     // for string
#include <vector>                     // for vector
#include "input/accessors.hh"         // for Record
#include "mesh/point.hh"
#include "tools/general_iterator.hh"  // for Iter
#include "system/armor.hh"            // for Armor::array
#include "system/index_types.hh"      // for LongIdx

class Mesh;
class Distribution;
class OutputElement;
class OutputMesh;
class OutputMeshDiscontinuous;
namespace Input { namespace Type { class Record; } }
template<class T> class ElementDataCache;
template<int> class ElementAccessor;


typedef Iter<OutputElement> OutputElementIterator;


/**
 * @brief Base class for Output mesh.
 *
 * Defines common members for Output mesh classes:
 *  - OutputMesh represents output mesh with continuous elements
 *  - OutputMeshDiscontinuous represents output mesh with discontinuous elements
 *
 * Making of output meshes and calling of their initialization methods must be execute in correct order, see example:
@code
    // Create or get Mesh object
    Mesh * my_mesh = ...

    // Construct mesh with continuous elements
    std::make_shared<OutputMesh> output_mesh = std::make_shared<OutputMesh>(*my_mesh);
    // Creates the sub meshes on all processes identical to the computational one.
    output_mesh->create_sub_mesh();
    // Creates mesh on zero process identical to the computational one.
    std::make_shared<OutputMesh> serial_output_mesh = output_mesh->make_serial_master_mesh();

    // Construct mesh with discontinuous elements
    std::make_shared<OutputMeshDiscontinuous> output_mesh_discont = std::make_shared<OutputMeshDiscontinuous>(*my_mesh);
    // Creates sub meshes on all processes mesh from the original my_mesh.
    output_mesh_discont->create_sub_mesh();
    // Creates mesh on zero process from the original my_mesh.
    std::make_shared<OutputMeshDiscontinuous> serial_output_mesh = output_mesh->make_serial_master_mesh();
@endcode
 */
class OutputMeshBase : public std::enable_shared_from_this<OutputMeshBase>
{
public:
    /// Shortcut instead of spacedim template. We suppose only spacedim=3 at the moment.
    static const unsigned int spacedim = 3;

    typedef std::function<void(const Armor::array &, const ElementAccessor<spacedim> &, std::vector<double> &)>
        ErrorControlFieldFunc;
    
    /// Constructor. Takes computational mesh as a parameter.
    OutputMeshBase(Mesh &mesh);
    /// Constructor. Takes computational mesh and input record as a parameters.
    OutputMeshBase(Mesh &mesh, const Input::Record &in_rec);
    virtual ~OutputMeshBase();
    
    /**
     * @brief The specification of output mesh.
     * @return record for output mesh
     */
    static const Input::Type::Record & get_input_type();
    
    /// Gives iterator to the FIRST element of the output mesh.
    OutputElementIterator begin();
    /// Gives iterator to the LAST element of the output mesh.
    OutputElementIterator end();
    
    /**
     * Creates sub mesh containing only local part of original (computation) mesh.
     *
     * TODO: should be replaced by local part of future parallel computational mesh.
     */
    void create_sub_mesh();

    /**
     * Creates refined sub mesh containing only local part of original (computation) mesh.
     */
    virtual void create_refined_sub_mesh()=0;

    /// Selects the error control field computing function of output field set according to input record.
    void set_error_control_field(ErrorControlFieldFunc error_control_field_func);

    /// Returns number of nodes.
    unsigned int n_nodes();
    /// Returns number of element.
    unsigned int n_elements();
    
    /// Check if nodes_, connectivity_ and offsets_ data caches are created
    bool is_created();

	/// Create nodes and elements data caches
	void create_id_caches();

	/// Synchronize parallel data and create serial COLECTIVE output mesh on zero process.
	void make_serial_master_mesh();

	/// Create output mesh of parallel output (implemented only for discontinuous mesh)
	virtual void make_parallel_master_mesh()
	{};

	/// Return master output mesh if exists or shared_ptr of this object.
	inline std::shared_ptr<OutputMeshBase> get_master_mesh() {
		if (master_mesh_) return master_mesh_;
		else return shared_from_this();
	};

protected:
	/**
	 * Possible types of OutputMesh.
	 */
	enum MeshType
	{
		orig,     //!< same as original (computational) mesh
		refined,  //!< refined mesh
		discont   //!< discontinuous mesh
	};

	/**
	 * Construct empty output mesh.
	 *
	 * Use in make_serial_master_mesh method and create mesh of same type as this object (continuous / discontinuos)
	 */
	virtual std::shared_ptr<OutputMeshBase> construct_mesh()=0;

	/**
	 * Create serial (collective) nodes cache on zero process.
	 *
	 * Implements part of \p make_serial_master_mesh that are specific for continuous and discontinuous case.
	 */
	virtual std::shared_ptr<ElementDataCache<double>> make_serial_nodes_cache(std::shared_ptr<ElementDataCache<unsigned int>> global_offsets)=0;

	/**
	 * Create serial (collective) connectivity cache on zero process.
	 *
	 * Implements part of \p make_serial_master_mesh that are specific for continuous and discontinuous case.
	 */
	virtual std::shared_ptr<ElementDataCache<unsigned int>> make_serial_connectivity_cache(std::shared_ptr<ElementDataCache<unsigned int>> global_offsets)=0;

	/// Compute and return number of nodes for each elements (compute from offsets)
	std::shared_ptr<ElementDataCache<unsigned int>> get_elems_n_nodes();

	/// Input record for output mesh.
    Input::Record input_record_;
    
    /// Pointer to the computational mesh.
    Mesh *orig_mesh_;
    
    /// Maximal level of refinement.
    const unsigned int max_level_;
    
    /// Refinement error control field function (hold value_list function of field).
    ErrorControlFieldFunc error_control_field_func_;

    MeshType mesh_type_;                ///< Type of OutputMesh
    bool refine_by_error_;              ///< True, if output mesh is to be refined by error criterion.
    double refinement_error_tolerance_; ///< Tolerance for error criterion refinement.
    
    /// Vector of element indices in the computational mesh. (Important when refining.)
    std::shared_ptr<std::vector<unsigned int>> orig_element_indices_;

    /// Vector of node coordinates. [spacedim x n_nodes]
    std::shared_ptr<ElementDataCache<double>> nodes_;
    /// Vector maps the nodes to their coordinates in vector @p nodes_.
    std::shared_ptr<ElementDataCache<unsigned int>> connectivity_;
    /// Vector of offsets of node indices of elements. Maps elements to their nodes in connectivity_.
    std::shared_ptr<ElementDataCache<unsigned int>> offsets_;

    /// Vector gets ids of nodes. Data is used in GMSH output.
    std::shared_ptr<ElementDataCache<unsigned int>> node_ids_;
    /// Vector gets ids of elements. Data is used in GMSH output.
    std::shared_ptr<ElementDataCache<unsigned int>> elem_ids_;
    /// Vector gets ids of regions. Data is used in GMSH output.
    std::shared_ptr<ElementDataCache<unsigned int>> region_ids_;
    /// Vector gets partitions of elements. Data is used in GMSH output.
    std::shared_ptr<ElementDataCache<int>> partitions_;

    /**
     * Master OutputMesh.
     *
     *  - serial output: is constructed on zero process (collective) and allow to produce serial output of parallel computation
     *  - parallel output: is constructed on each process only for discontinuous mesh
     */
    std::shared_ptr<OutputMeshBase> master_mesh_;

    /**
     * Next variables hold distributions of elements and nodes. They differ for mesh types
     *  - continuous and discontinuous mesh shared objects with computational (orig) mesh
     *  - refined mesh creates own objects
     */
    LongIdx *el_4_loc_;           ///< Index set assigning to local element index its global index.
    Distribution *el_ds_;         ///< Parallel distribution of elements.
    LongIdx *node_4_loc_;         ///< Index set assigning to local node index its global index.
    Distribution *node_ds_;       ///< Parallel distribution of nodes. Depends on elements distribution.
    unsigned int n_local_nodes_;  ///< Hold number of local nodes (own + ghost), value is equal with size of node_4_loc array.

    /// Friend provides access to vectors for element accessor class.
    friend class OutputElement;
    friend class OutputTime;
    friend class OutputMSH;
    friend class OutputVTK;
    friend class OutputMesh;
    friend class OutputMeshDiscontinuous;
};


/// @brief Class represents output mesh with continuous elements.
class OutputMesh : public OutputMeshBase
{
public:
    OutputMesh(Mesh &mesh);
    OutputMesh(Mesh &mesh, const Input::Record &in_rec);
    ~OutputMesh();
    
    /// Implements OutputMeshBase::create_refined_sub_mesh
    void create_refined_sub_mesh() override;

protected:
    bool refinement_criterion();
    
    /// Implements OutputMeshBase::construct_mesh
    std::shared_ptr<OutputMeshBase> construct_mesh() override;

    /// Implements OutputMeshBase::make_serial_nodes_cache
    std::shared_ptr<ElementDataCache<double>> make_serial_nodes_cache(std::shared_ptr<ElementDataCache<unsigned int>> global_offsets) override;

    /// Implements OutputMeshBase::make_serial_connectivity_cache
	std::shared_ptr<ElementDataCache<unsigned int>> make_serial_connectivity_cache(std::shared_ptr<ElementDataCache<unsigned int>> global_offsets) override;

	/// Friend provides access to vectors for discontinous output mesh.
    friend class OutputMeshDiscontinuous;
};


/// @brief Class represents output mesh with discontinuous elements.
class OutputMeshDiscontinuous : public OutputMeshBase
{
public:
    OutputMeshDiscontinuous(Mesh &mesh);
    OutputMeshDiscontinuous(Mesh &mesh, const Input::Record& in_rec);
    ~OutputMeshDiscontinuous();
    
    /// Implements OutputMeshBase::create_refined_sub_mesh
    void create_refined_sub_mesh() override;

    /// Overrides OutputMeshBase::make_parallel_master_mesh
    void make_parallel_master_mesh() override;

protected:
    ///Auxiliary structure defining element of refined output mesh.
    struct AuxElement{
        std::vector<Space<spacedim>::Point> nodes;
        unsigned int level;
    };

    ///Performs the actual refinement of AuxElement. Recurrent.
    template<int dim>
    void refine_aux_element(const AuxElement& aux_element,
                            std::vector< AuxElement >& refinement,
                            const ElementAccessor<spacedim> &ele_acc
                           );
    
    /// Collects different refinement criteria results.
    bool refinement_criterion(const AuxElement& ele,
                              const ElementAccessor<spacedim> &ele_acc);
    
    /// Refinement flag - checks only maximal level of refinement.
    bool refinement_criterion_uniform(const AuxElement& ele);
    
    /// Refinement flag - measures discretisation error according to error control field.
    bool refinement_criterion_error(const AuxElement& ele,
                                    const Space<spacedim>::Point &centre,
                                    const ElementAccessor<spacedim> &ele_acc
                                   );

    /// Implements OutputMeshBase::construct_mesh
    std::shared_ptr<OutputMeshBase> construct_mesh() override;

    /// Implements OutputMeshBase::make_serial_nodes_cache
    std::shared_ptr<ElementDataCache<double>> make_serial_nodes_cache(std::shared_ptr<ElementDataCache<unsigned int>> global_offsets) override;

    /// Implements OutputMeshBase::make_serial_connectivity_cache
	std::shared_ptr<ElementDataCache<unsigned int>> make_serial_connectivity_cache(std::shared_ptr<ElementDataCache<unsigned int>> global_offsets) override;
};

#endif  // OUTPUT_MESH_HH_

