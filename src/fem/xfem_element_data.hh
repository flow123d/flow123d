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
#include "quadrature/qxfem.hh"
#include "fem/global_enrichment_func.hh"

template<int dim, int spacedim> class QXFEM;

enum Quantity{
    velocity = 0,
    pressure = 1,
    pressure_lagrange = 2
};

/** @brief Base class for XFEM element data.
 */
class XFEMElementDataBase
{
public:
    /// Constructor.
    XFEMElementDataBase()
    : ele_global_idx_(-1)
    {}
    
    ///Destructor
    virtual ~XFEMElementDataBase()
    {}
    
    ///Index of enriched element
    unsigned int ele_global_idx() const
    { return ele_global_idx_;}
    
    ///Index of element intersecting this enriched element
    unsigned int intersection_ele_global_idx() const
    { return intersection_ele_global_idx_;}
    
    void set_element(int global_idx, int intersection_global_idx)
    { ele_global_idx_ = global_idx;
      intersection_ele_global_idx_ = intersection_global_idx;
    }
    
    ///Returns number of wells comunicating with the cell
    unsigned int n_enrichments() const
    { return global_enrichment_indices_.size();}
    
    /// Returns global index of the well.
    /** @param local_well_index is local well index in the cell
     */
    unsigned int global_enrichment_index(unsigned int local_enrichment_index) const
    {   ASSERT_DBG(local_enrichment_index < n_enrichments());
        return global_enrichment_indices_[local_enrichment_index];}
        
    const std::vector<unsigned int> & global_enrichment_indices() const
    {   return global_enrichment_indices_;}
    
    
    /// Getter for enriched dofs vectors: [quantity][local_enrichment_index][local_dof].
    std::vector<std::vector<std::vector<int>>> &global_enriched_dofs()
    { return global_enriched_dofs_;}
    
    /// Getter for enriched dofs by a single quantity and a single well.
    const std::vector<int> &global_enriched_dofs(Quantity quant,
                                                 unsigned int local_enrichment_index) const;
    
    /// Number of all degrees of freedom on the cell (from all quantities, all enrichments).
    unsigned int n_enriched_dofs() const;
      
    /// Number of degrees of freedom from a single quantity @p quant, all enrichments.
    unsigned int n_enriched_dofs(Quantity quant) const;
      
    /// Number of degrees of freedom from a single @p quant, a single enrichment @p local_enrichment_index.
    unsigned int n_enriched_dofs(Quantity quant, unsigned int local_enrichment_index) const;
    
    /** Returns true if @p local_enrichment_index enrichment has nonzero cross-section with the element.
     */
    bool enrichment_intersects(unsigned int local_enrichment_index) const
    {   ASSERT_DBG(local_enrichment_index < n_enrichments());
        return enrichment_intersects_[local_enrichment_index];}
        
    /** Number of wells that has nonzero cross-section with the element.
     * It is determined by non-zero count of quadrature points along the singularity edge.
     */
    unsigned int n_enrichments_intersect() const;
    
    void print(std::ostream& out) const;

protected:
///iterator of the cell to which this data object belongs
//     LocalElementAccessorBase<spacedim> ele_ac_;
//     ElementFullIter ele_;
    
    /// Index of element to which this data belongs to.
    unsigned int ele_global_idx_;
    
    unsigned int intersection_ele_global_idx_;
    
    ///global indices of the wells
    std::vector<unsigned int> global_enrichment_indices_;
    
    ///true for enrichments that intersects the element
    std::vector<bool> enrichment_intersects_;
    
    /** Global numbers of enriched DoFs. 
     * Index subset in \f$ \mathcal{M}_w \f$ (nodes on both reproducing and blending elements).
     * Access the index in format [quantity][enrichment_index][local_node_index].
     */
    std::vector<std::vector<std::vector<int>>> global_enriched_dofs_;
};


/** @brief Base class for data distributed umong enriched elements.
 * These are pointers to @p GlobalEnrichmentFunc objects, quadrature points of wells,
 * and enriched degrees of freedom.
 */
template<int dim, int spacedim>
class XFEMElementData : public XFEMElementDataBase
{
public:
    typedef typename std::shared_ptr<GlobalEnrichmentFunc<dim,spacedim>> EnrichmentPtr;
    
    /// Constructor.
    XFEMElementData()
    {}
    
    ///Destructor
    virtual ~XFEMElementData()
    {}

    /// @name Getters
    //@{
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
    
    
//     /** Writes local DoFs in given vector: wells*[FE dofs, Xdofs, Wdofs]
//      * Sets n_wells_inside, n_dofs, n_xdofs, n_wdofs.
//      */
//     void get_dof_indices(std::vector<int> &local_dof_indices, unsigned int fe_dofs_per_cell);
    
    
    /** Getter for enrichment function value of a single well at nodes.
     * Provides acces to the map of node values of enrichment functions.
     */
//     double node_enrich_value(unsigned int local_enrichment_index, unsigned int local_vertex_index) const;
    //@}
    
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
    
//     void set_node_values(std::vector<std::map<int, double> > *node_vals,
//                          std::vector<std::map<int, Space<3>::Point> > *node_vec_vals)
//     {   node_values = node_vals; 
//         node_vec_values = node_vec_vals;
//     }
  
protected:   
    ///vector of pointers to wells
    std::vector<EnrichmentPtr> enrichment_func_;
    
//     std::vector<std::map<int, double> > *node_values;
//     std::vector<std::map<int, Space<3>::Point> > *node_vec_values;
};
    

#include "fem/singularity.hh"
//*************************************************************************************
//*************************************************************************************

/** @brief Class storing data from wells distributed to cells. Used in class @p XModel.
 * 
 * This class stores similar data as class @p DataCell
 * but also the enriched degrees of freedom belonging to the current cell.
 */
template<int dim>
class XFEMElementSingularData : public XFEMElementData<dim,3>
{
  public:
    
    /// Constructor.
    XFEMElementSingularData();
        
    /// Destructor
    virtual ~XFEMElementSingularData();
    
    const QXFEM<dim,3>& sing_quadrature(unsigned int local_enrichment_index) const;
    
    void create_sing_quads(ElementFullIter ele);
    
    //TODO: get rid of this hack
    std::vector<std::shared_ptr<Singularity<dim-2>>> sing_vec(){
        std::vector<std::shared_ptr<Singularity<dim-2>>> res(this->enrichment_func_.size());
        for(unsigned int i=0; i<this->enrichment_func_.size(); i++){
            res[i] = std::static_pointer_cast<Singularity<dim-2>>(this->enrichment_func_[i]);
        }
        return res;
    }
  private:
    
    std::vector<QXFEM<dim,3>> sing_quads_;
};


#endif // XFEM_ELEMENT_DATA_HH_
