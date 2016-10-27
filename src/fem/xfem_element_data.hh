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
 * @file    xfem_element_data.hh
 * @brief   XFEM data class associated with an enriched element.
 * @author  Pavel Exner
 */
#ifndef XFEM_ELEMENT_DATA_HH_
#define XFEM_ELEMENT_DATA_HH_

#include "mesh/mesh.h"
// #include "flow/mh_dofhandler.hh"
#include "fem/global_enrichment_func.hh"
#include "quadrature/qxfem.hh"

class Well;
template<int dim, int spacedim> class QXFEM;

/** @brief Base class for data distributed umong enriched elements.
 * These are pointers to @p GlobalEnrichmentFunc objects, quadrature points of wells,
 * and enriched degrees of freedom.
 */
template<int dim, int spacedim>
class XFEMElementDataBase
{
public:
    typedef typename std::shared_ptr<GlobalEnrichmentFunc<dim,spacedim>> EnrichmentPtr;
    
    /// Constructor.
    XFEMElementDataBase()
    {}
    
    ///Destructor
    virtual ~XFEMElementDataBase()
    {}

    /// @name Getters
    //@{
    unsigned int ele_global_idx() const
    { return ele_global_idx_;}
    
    void set_element(int global_idx)
    { ele_global_idx_ = global_idx;}
    
    ///Returns number of wells comunicating with the cell
    unsigned int n_enrichments() const
    { return enrichment_func_.size();}
    
    /// Returns pointer to one of the wells comunicating with the cell this data belong to.
    /**
     * @param local_well_index is local well index in the cell
     * @return pointer to well
     */
    EnrichmentPtr enrichment_func(unsigned int local_enrichment_index) const
    {   ASSERT_DBG(local_enrichment_index < n_enrichments());
        return enrichment_func_[local_enrichment_index];}
    
     /// Returns pointer to one of the wells comunicating with the cell this data belong to.
    /**
     * @param local_well_index is local well index in the cell
     * @return constant reference to a vector of pointers to wells
     */
    const std::vector<EnrichmentPtr> & enrichment_func_vec() const
    { return enrichment_func_;}
    
    /// Returns global index of the well.
    /** @param local_well_index is local well index in the cell
     */
    unsigned int global_enrichment_index(unsigned int local_enrichment_index) const
    {   ASSERT_DBG(local_enrichment_index < n_enrichments());
        return global_enrichment_indices_[local_enrichment_index];}
        
    /// Returns global dof index of the well.
    /** @param local_well_index is local well index in the cell_
     */
//     unsigned int get_well_dof_index(const unsigned int &local_well_index);
    
    ///Returns reference to vector of pointers to quadrature points of the well boundary
//     const std::vector<const dealii::Point<2>* > &q_points(unsigned int local_well_index);
    
    ///Returns reference to vector of pointers to quadrature points of the well boundary
//     const std::vector<dealii::Point<2> > &mapped_q_points(unsigned int local_well_index);
    //@}
    
    /// Maps the quadrature points lying in the cell to a reference cell.
//     void map_well_quadrature_points(const dealii::Mapping<2>& mapping);
    
    /// Sets global well indices of the wells in the cell.
    /** @param well_dof_indices is vector of global well indices.
     */
//     void set_well_dof_indices(const std::vector<unsigned int> &well_dof_indices);
    
    /// Adds new data to this object.
    /**
     * @param well is pointer to well which lies in the cell
     * @param well_index is index of the well in the global vector of wells in model class
     */
    void add_data(EnrichmentPtr enrichment_func,
                  unsigned int global_enrichment_index)
    {
        enrichment_func_.push_back(enrichment_func);
        global_enrichment_indices_.push_back(global_enrichment_index);
    }
    
    void set_node_values(std::vector<std::map<int, double> > *node_vals,
                         std::vector<std::map<int, Space<3>::Point> > *node_vec_vals)
    {   node_values = node_vals; 
        node_vec_values = node_vec_vals;
    }
  
protected:
    ///iterator of the cell to which this data object belongs
//     LocalElementAccessorBase<spacedim> ele_ac_;
//     ElementFullIter ele_;
    unsigned int ele_global_idx_;
    
    ///vector of pointers to wells
    std::vector<EnrichmentPtr> enrichment_func_;
    ///global indices of the wells
    std::vector<unsigned int> global_enrichment_indices_;
    ///global dof indices of the wells
//     std::vector<unsigned int> well_dof_indices_;
    
    std::vector<std::map<int, double> > *node_values;
    
    std::vector<std::map<int, Space<3>::Point> > *node_vec_values;
  
    /** Pointers to quadrature points of wells that lies inside the cell.
     * Access: Point = [local_well_index][q] 
     */
//     std::vector< std::vector< const dealii::Point<2>* > > q_points_;
    
    /// Mapped well quadrature points (@p q_points_ on a reference cell).
//     std::vector< std::vector<dealii::Point<2> > > mapped_q_points_;
    
    ///just for returning zero lenght vector
//     std::vector< const dealii::Point<2>* > dummy_q_points_;
};



enum Quantity{
    velocity = 0,
    pressure = 1,
    pressure_lagrange = 2
};
    
//*************************************************************************************
//*************************************************************************************

/** @brief Class storing data from wells distributed to cells. Used in class @p XModel.
 * 
 * This class stores similar data as class @p DataCell
 * but also the enriched degrees of freedom belonging to the current cell.
 */
// template<int dim, int spacedim>
class XFEMElementSingularData : public XFEMElementDataBase<2,3>
{
  public:
     
    
    /// Constructor.
    XFEMElementSingularData();
    
//     XFEMElementSingularData(LocalElementAccessorBase<3> ele_ac)
//       : XFEMElementDataBase<2,3>(ele_ac)
//     XFEMElementSingularData(ElementFullIter ele)
//     : XFEMElementDataBase<2,3>(ele)
//     {}
    
//     /// Constructor. 
//     /// For using without quadrature points around the well edge.
//     XFEMElementData(LocalElementAccessorBase<3> ele_ac, 
//               Well *well, 
//               const unsigned int &well_index,
//               const std::vector<unsigned int> &enriched_dofs,
//               const std::vector<unsigned int> &weights);
//     
//     /// Constructor. 
//     /// For using with quadrature points around the well edge.
//     XFEMElementData(LocalElementAccessorBase<3> ele_ac, 
//               Well *well, 
//               const unsigned int &well_index,
//               const std::vector< unsigned int > &enriched_dofs,
//               const std::vector<unsigned int> &weights,
//               const std::vector<const dealii::Point<2>* > &q_points);
    
    /// Destructor
    virtual ~XFEMElementSingularData();
    
    /// @name Getters
    //@{    
    /// Getter for enriched dofs vectors: [quantity][local_enrichment_index][local_dof].
    std::vector<std::vector<std::vector<int>>> &global_enriched_dofs();
    
     
    /// Getter for enriched dofs by a single quantity and a single well.
    const std::vector<int> &global_enriched_dofs(Quantity quant,
                                                 unsigned int local_enrichment_index) const;
       
      /// Getter for weights of a single well.
//       const std::vector<unsigned int> &weights(unsigned int local_enrichment_index);
      
      
      /** Getter for enrichment function value of a single well at nodes.
       * Provides acces to the map of node values of enrichment functions.
       */
      double node_enrich_value(unsigned int local_enrichment_index, unsigned int local_vertex_index) const;
      
      /** Writes local DoFs in given vector: wells*[FE dofs, Xdofs, Wdofs]
       * Sets n_wells_inside, n_dofs, n_xdofs, n_wdofs.
       */
      void get_dof_indices(std::vector<int> &local_dof_indices, unsigned int fe_dofs_per_cell);
      
      /// Number of all degrees of freedom on the cell (from all quantities, all enrichments).
      unsigned int n_enriched_dofs() const;
      
    /// Number of degrees of freedom from a single quantity @p quant, all enrichments.
    unsigned int n_enriched_dofs(Quantity quant) const;
      
    /// Number of degrees of freedom from a single @p quant, a single enrichment @p local_enrichment_index.
    unsigned int n_enriched_dofs(Quantity quant, unsigned int local_enrichment_index) const;
      
    /** Number of wells that has nonzero cross-section with the element.
     * It is determined by non-zero count of quadrature points along the singularity edge.
     */
    unsigned int n_singularities_inside() const;
    
    /** Returns true if @p local_enrichment_index singularity has nonzero cross-section with the element.
     * It is determined by non-zero count of quadrature points along the singularity edge.
     */
    bool is_singularity_inside(unsigned int local_enrichment_index) const;
      
      /// Number of all degrees of freedom on the cell.
//       unsigned int n_standard_dofs();
      
      /// Number of all degrees of freedom on the cell.
//       unsigned int n_dofs();
      
    const QXFEM<2,3>& sing_quadrature(unsigned int local_enrichment_index) const;
      
      /// Number of polar quadratures for wells.
//       unsigned int n_polar_quadratures(void);
      
//       XQuadratureWell * polar_quadrature(unsigned int local_well_index);
//       std::vector<XQuadratureWell *> polar_quadratures(void);
    //@}
    
    /// Add enriched data (without q_points).
//     void add_data(Well *well, 
//                   const unsigned int &well_index, 
//                   const std::vector<unsigned int> &enriched_dofs,
//                   const std::vector<unsigned int> &weights);
    
    /// Add enriched data (possibly with q_points).
//     void add_data(Well *well, 
//                   const unsigned int &well_index, 
//                   const std::vector<unsigned int> &enriched_dofs,
//                   const std::vector<unsigned int> &weights,
//                   const std::vector<const dealii::Point<2>* > &q_points);
    
//     void set_polar_quadrature(XQuadratureWell* xquad);
    
//     void clear_polar_quadratures(void);
    
    void create_sing_quads(ElementFullIter ele);
    
//     /** STATIC function. Goes through given XFEMElementSingularDatas objects and initialize node values of enrichment before system assembly.
//      * @param data_vector is given output vector (by wells) of maps which map enrichment values to the nodes
//      * @param xdata is given vector of XFEMElementSingularData objects (includes enrichment functions and cells)
//      * @param n_wells is the total number of wells in the model
//      */
//     static void initialize_node_values(std::vector<std::map<unsigned int, double> > &data_vector, 
//                                        std::vector<XFEMElementSingularData*> xdata, 
//                                        unsigned int n_wells);
    
    void print(std::ostream& out) const;
  private:
    
    std::vector<QXFEM<2,3>> sing_quads_;
    
    /// Quadratures in polar coordinates in vicinity of wells affecting the current cell.
//     std::vector<XQuadratureWell*> well_xquadratures_;
    
    /** Global numbers of enriched DoFs. 
     * Index subset in \f$ \mathcal{M}_w \f$ (nodes on both reproducing and blending elements).
     * Access the index in format [quantity][well_index][local_node_index].
     */
    std::vector<std::vector<std::vector<int>>> global_enriched_dofs_;
    
//     std::vector<unsigned int> n_enriched_dofs_per_well_; ///<Number of enriched dofs by a single well.
//     unsigned int 
//                 n_enriched_dofs_,              ///< Number of all enriched dofs.
//                  n_wells_inside_,               ///< Number of wells inside the cell.
//                  n_standard_dofs_,              ///< Number of standard dofs.
//                  n_dofs_,                       ///< Total number of dofs.
//                  n_polar_quadratures_;          ///< Number of polar quadratures for wells.
    
    /** Weights of enriched nodes. 
     * Weight is equal \f$ g_u = 1 \$ at enriched node from subset \f$ \mathcal{N}_w \f$.
     * Weight is equal \f$ g_u = 0 \$ at enriched node from subset \f$ \mathcal{M}_w \f$ which is not in \f$ \mathcal{N}_w \f$
     * Access the index in format [local_well_index][local_node_index].
     */
//     std::vector<std::vector<unsigned int> > weights_;
};


#endif // XFEM_ELEMENT_DATA_HH_