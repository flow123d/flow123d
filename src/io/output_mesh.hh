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

#include <string>
#include <functional>

#include "system/sys_profiler.hh"
#include "input/accessors.hh"

#include "tools/general_iterator.hh"
#include "mesh/point.hh"

class Mesh;
template<class T> class ElementDataCache;
template<int> class ElementAccessor;

class OutputElement;
typedef GeneralIterator<OutputElement> OutputElementIterator;

class OutputMeshBase;
class OutputMesh;
class OutputMeshDiscontinuous;
class OutputMSH;
class OutputVTK;


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
    // Creates the mesh identical to the computational one.
    output_mesh->create_mesh();

    // Construct mesh with discontinuous elements
    std::make_shared<OutputMeshDiscontinuous> output_mesh_discont = std::make_shared<OutputMeshDiscontinuous>(*my_mesh);
    // Creates mesh from the original my_mesh.
    output_mesh_discont->create_mesh();
@endcode
 */
class OutputMeshBase : public std::enable_shared_from_this<OutputMeshBase>
{
public:
    /// Shortcut instead of spacedim template. We suppose only spacedim=3 at the moment.
    static const unsigned int spacedim = 3;

    typedef std::function<void(const std::vector< Space<spacedim>::Point > &, const ElementAccessor<spacedim> &, std::vector<double> &)>
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
    
    /// Creates the output mesh identical to the orig mesh.
    virtual void create_mesh()=0;

    /// Creates refined mesh.
    virtual void create_refined_mesh()=0;

    /// Creates sub mesh containing only local elements.
    virtual void create_sub_mesh()=0;

    /// Selects the error control field computing function of output field set according to input record.
    void set_error_control_field(ErrorControlFieldFunc error_control_field_func);

    /// Returns number of nodes.
    unsigned int n_nodes();
    /// Returns number of element.
    unsigned int n_elements();
    
    /// Check if nodes_, connectivity_ and offsets_ data caches are created
    bool is_created();

    /// Return data cache of node ids. If doesn't exist create its.
    std::shared_ptr<ElementDataCache<unsigned int>> get_node_ids_cache();
    /// Return data cache of element ids. If doesn't exist create its.
	std::shared_ptr<ElementDataCache<unsigned int>> get_element_ids_cache();

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


	/// Create node_ids_ and elem_ids_ data caches
	void create_id_caches();


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

    /// Friend provides access to vectors for element accessor class.
    friend class OutputElement;
    friend class OutputMSH;
    friend class OutputVTK;
};


/// @brief Class represents output mesh with continuous elements.
class OutputMesh : public OutputMeshBase
{
public:
    OutputMesh(Mesh &mesh);
    OutputMesh(Mesh &mesh, const Input::Record &in_rec);
    ~OutputMesh();
    
    /// Creates the output mesh identical to the orig mesh.
    void create_mesh() override;
    
    /// Creates refined mesh.
    void create_refined_mesh() override;
    
    /// Creates sub mesh.
    void create_sub_mesh() override;

protected:
    bool refinement_criterion();
    
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
    
    /// Creates the output mesh identical to the orig mesh.
    void create_mesh() override;
    
    /// Creates discontinuous refined mesh.
    void create_refined_mesh() override;
    
    /// Creates sub mesh.
    void create_sub_mesh() override;

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
};

#endif  // OUTPUT_MESH_HH_

