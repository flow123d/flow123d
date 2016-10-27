#ifndef XFEM_ELEMENT_DATA_HH_
#define XFEM_ELEMENT_DATA_HH_

#include "mesh/mesh.h"
#include "flow/mh_dofhandler.hh"
// #include <flow/darcy_flow_xfem.hh>
#include "fem/global_enrichment_func.hh"

class XQuadratureWell;

//forward declarations
// namespace dealii{
//     template<int,int> class Mapping;
// }

class Well;
template<int dim, int spacedim> class QXFEM;

/** @brief Base class for data distributed umong cells.
 * We need to distribute some data from wells umong the cells
 * These are pointers to @p Well objects, quadrature points of wells,
 * in case of XFEM additional enriched degrees of freedom.
 */
template<int dim, int spacedim>
class XFEMElementDataBase
{
public:
    typedef typename std::shared_ptr<GlobalEnrichmentFunc<dim,spacedim>> EnrichmentPtr;
    
    XFEMElementDataBase()
    {}
    
    /** @brief Constructor.
     * @param cell iterator to cell which this data belongs to
     */
//     XFEMElementDataBase(LocalElementAccessorBase<spacedim> ele_ac)
//       : ele_ac_(ele_ac)
//     XFEMElementDataBase(ElementFullIter ele)
//     : ele_(ele)
//     {
//         ASSERT_DBG(ele->dim() == dim);
//     }

    ///Destructor
    virtual ~XFEMElementDataBase()
    {}

    /// @name Getters
    //@{
    /// Returns pointer to the cell which this data belong to.
//     LocalElementAccessorBase<spacedim> element();
//     ElementFullIter element()
//     { return ele_;}
    unsigned int ele_global_idx()
    { return ele_global_idx_;}
    
    void set_element(int global_idx)
    { ele_global_idx_ = global_idx;}
    
    ///Returns number of wells comunicating with the cell
    unsigned int n_enrichments()
    { return enrichment_func_.size();}
    
    /// Returns pointer to one of the wells comunicating with the cell this data belong to.
    /**
     * @param local_well_index is local well index in the cell
     * @return pointer to well
     */
    EnrichmentPtr enrichment_func(unsigned int local_enrichment_index)
    {   ASSERT_DBG(local_enrichment_index < n_enrichments());
        return enrichment_func_[local_enrichment_index];}
    
     /// Returns pointer to one of the wells comunicating with the cell this data belong to.
    /**
     * @param local_well_index is local well index in the cell
     * @return constant reference to a vector of pointers to wells
     */
    const std::vector<EnrichmentPtr> & enrichment_func_vec()
    { return enrichment_func_;}
    
    /// Returns global index of the well.
    /** @param local_well_index is local well index in the cell
     */
    unsigned int global_enrichment_index(unsigned int local_enrichment_index)
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
    void add_data(EnrichmentPtr enrichment_func, unsigned int global_enrichment_index)
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
    
//     unsigned int n_vertices_;             ///< Number of vertices.
  
    /** Pointers to quadrature points of wells that lies inside the cell.
     * Access: Point = [local_well_index][q] 
     */
//     std::vector< std::vector< const dealii::Point<2>* > > q_points_;
    
    /// Mapped well quadrature points (@p q_points_ on a reference cell).
//     std::vector< std::vector<dealii::Point<2> > > mapped_q_points_;
    
    ///just for returning zero lenght vector
//     std::vector< const dealii::Point<2>* > dummy_q_points_;
};


//*************************************************************************************
//*************************************************************************************

/** @brief Class storing data from wells distributed to cells. Used in class @p Model.
 * 
 * This class is used to store data at the cell in the class @p Model.
 */
// class DataCell : public XFEMElementDataBase
// {
// public:
//   ///Constructor.
//   /**
//    * @param cell iterator to cell which this data belongs to
//    */
//   DataCell(LocalElementAccessorBase<3> element) 
//     : XFEMElementDataBase(element)
//   {}
//   
//   ///Constructor.
//   /**
//    * @param cell iterator to cell which this data belongs to
//    * @param well is pointer to well which lies in the cell
//    * @param well_index is index of the well in the global vector of wells in model class
//    * @param q_points is vector of pointers to quadrature points of the well that lie in the cell
//    * */
//   DataCell(LocalElementAccessorBase<spacedim> ele_ac, 
//             Well *well, 
//             const unsigned int &well_index,
//             const std::vector<const dealii::Point<2>* > &q_points);
//   
//   ///Destructor.
//   virtual ~DataCell()
//   {}
//       
//   /// Adds new data to this object.
//   /**
//    * @param well is pointer to well which lies in the cell
//    * @param well_index is index of the well in the global vector of wells in model class
//    * @param q_points is vector of pointers to quadrature points
//    */
//   void add_data(Well* well, 
//                 const unsigned int &well_index, 
//                 const std::vector<const dealii::Point<2>* > &q_points);
// };


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
    XFEMElementSingularData()
    {}
    
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
    virtual ~XFEMElementSingularData()
    {}
    
    /// @name Getters
    //@{    
    /// Getter for enriched dofs by a single well.
    std::vector<std::vector<std::vector<int>>> &global_enriched_dofs()
    { return global_enriched_dofs_;}
     
    /// Getter for enriched dofs by a single well.
    const std::vector<int> &global_enriched_dofs(Quantity quant,
                                                 unsigned int local_enrichment_index)
    {   ASSERT_DBG(quant < global_enriched_dofs_.size());
        ASSERT_DBG(local_enrichment_index < global_enriched_dofs_[quant].size());
        return global_enriched_dofs_[quant][local_enrichment_index];
    }
       
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
      unsigned int n_enriched_dofs();
      
    /// Number of degrees of freedom on the cell (from a single quantity @p quant, all enrichments).
    unsigned int n_enriched_dofs(Quantity quant){
        ASSERT_DBG(quant < global_enriched_dofs_.size());
        ASSERT_DBG(global_enriched_dofs_[quant].size() > 0);
        ASSERT_DBG(global_enriched_dofs_[quant][0].size() > 0);
        //FIXME: supposing that all enrichments have the same number of enr dofs
        return global_enriched_dofs_[quant].size() * global_enriched_dofs_[quant][0].size();
    }
      
    /// Number of degrees of freedom on the cell (from a single @p quant, a single enrichment @p local_enrichment_index).
    unsigned int n_enriched_dofs(Quantity quant, unsigned int local_enrichment_index){
        ASSERT_DBG(quant < global_enriched_dofs_.size());
        ASSERT_DBG(local_enrichment_index < global_enriched_dofs_[quant].size());
        return global_enriched_dofs_[quant][local_enrichment_index].size();
    }
      
      /// Number of wells that has nonzero cross-section with the cell.
      unsigned int n_wells_inside();
      
      /// Number of all degrees of freedom on the cell.
      unsigned int n_standard_dofs();
      
      /// Number of all degrees of freedom on the cell.
      unsigned int n_dofs();
      
    const QXFEM<2,3>& sing_quadrature(unsigned int local_enrichment_index)
    {   ASSERT_DBG(local_enrichment_index < sing_quads_.size());
        return sing_quads_[local_enrichment_index];
    }
      
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
    
    void clear_polar_quadratures(void);
    
    void create_sing_quads(ElementFullIter ele);
    
    /** STATIC function. Goes through given XFEMElementSingularDatas objects and initialize node values of enrichment before system assembly.
     * @param data_vector is given output vector (by wells) of maps which map enrichment values to the nodes
     * @param xdata is given vector of XFEMElementSingularData objects (includes enrichment functions and cells)
     * @param n_wells is the total number of wells in the model
     */
    static void initialize_node_values(std::vector<std::map<unsigned int, double> > &data_vector, 
                                       std::vector<XFEMElementSingularData*> xdata, 
                                       unsigned int n_wells);
    
    void print(std::ostream& out);
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


#include "fem/singularity.hh"
#include "quadrature/qxfem.hh"
#include <armadillo>

inline void XFEMElementSingularData::create_sing_quads(ElementFullIter ele)
{
    const unsigned int n_qpoints = 100;
//     ElementFullIter ele = mesh_->element(xdata.ele_global_idx());
    sing_quads_.resize(n_enrichments());
    
    arma::mat proj = ele->element_map();
    
    std::map<unsigned int, arma::vec> unit_points_inside;
    
    DBGCOUT(<< "create_sing_quads on ele " << ele->index() << "\n");
    for(unsigned int w=0; w < n_enrichments(); w++){
        std::shared_ptr<Singularity0D<3>> sing = static_pointer_cast<Singularity0D<3>>(enrichment_func(w));
        sing->evaluate_q_points(n_qpoints);
        
        unit_points_inside.clear();
        
//         DBGCOUT(<< "test q_points\n");
        for(unsigned int q=0; q < n_qpoints; q++){
            const Space<3>::Point & p = sing->q_points()[q];
            arma::vec unit_p = ele->project_point(p, proj);
            
//             if(ele->index() == 42){
//                 sing->q_points()[q].print(cout,"real_p");
//                 unit_p.print(cout,"unit_p");
//             }
            
            if( unit_p(0) >= 0 && unit_p(0) <= 1 &&
                unit_p(1) >= 0 && unit_p(1) <= 1 &&
                unit_p(2) >= 0 && unit_p(2) <= 1){
        
//                 DBGCOUT(<< "qpoint inside\n");
                unit_points_inside[q] = unit_p;
            }
        }
        
        QXFEM<2,3>& qxfem = sing_quads_[w];
        qxfem.resize(unit_points_inside.size());
        std::map<unsigned int, arma::vec>::const_iterator pair;
        for(pair = unit_points_inside.begin(); pair != unit_points_inside.end(); pair++){
            
            qxfem.set_point(pair->first, RefElement<2>::bary_to_local(pair->second));
            qxfem.set_real_point(pair->first, sing->q_points()[pair->first]);
        }
        DBGCOUT(<< "quad[" << global_enrichment_index(w) << "] size " << sing_quads_[w].size() << "\n");
    }
}


inline void XFEMElementSingularData::print(ostream& out)
{
    out << this << "xdata: ele " << ele_global_idx_ << " enrichments " << n_enrichments();
    out << " dofs[ ";
//     for(unsigned int q=0; q<global_enriched_dofs_.size(); q++)
    for(unsigned int q=0; q<2; q++)
        for(unsigned int w=0; w<n_enrichments(); w++)
            for(unsigned int j=0; j<global_enriched_dofs_[q][w].size(); j++){
                out << global_enriched_dofs_[q][w][j] << " ";
            }
    out << "]\n";
}






/****************************************            Implementation          ********************************/
/*
template<int dim, int spacedim>
// inline LocalElementAccessorBase<spacedim> XFEMElementDataBase<dim,spacedim>::element() { return ele_ac_; }
inline ElementFullIter XFEMElementDataBase<dim,spacedim>::element() { return ele_; }

template<int dim, int spacedim>
inline unsigned int XFEMElementDataBase<dim,spacedim>::n_enrichments() { return enrichment_func_.size(); } */

// template<int dim, int spacedim>
// inline const std::vector< GlobalEnrichmentFunc<dim,spacedim>* >& XFEMElementDataBase<dim,spacedim>::get_wells() { return enrichment_func_;}

#endif // XFEM_ELEMENT_DATA_HH_