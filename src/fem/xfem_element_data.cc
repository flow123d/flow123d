#include "xfem_element_data.hh"
#include "singularity.hh"
#include "quadrature/qxfem.hh"
#include "quadrature/qxfem_factory.hh"

#include <armadillo>
#include "mapping_p1.hh"

using namespace std;

const vector< int >& XFEMElementDataBase::global_enriched_dofs(Quantity quant,
                                                               unsigned int local_enrichment_index) const
{   ASSERT_DBG(quant < global_enriched_dofs_.size());
    ASSERT_DBG(local_enrichment_index < global_enriched_dofs_[quant].size());
    return global_enriched_dofs_[quant][local_enrichment_index];
}

unsigned int XFEMElementDataBase::n_enriched_dofs(Quantity quant) const
{
//     DBGVAR(ele_global_idx_); cout<<endl;
//     DBGVAR(global_enriched_dofs_.size()); cout<<endl;
//     DBGVAR(global_enriched_dofs_[quant].size()); cout<<endl;
//     DBGVAR(global_enriched_dofs_[quant][0].size()); cout<<endl;
    
    ASSERT_DBG(quant < global_enriched_dofs_.size());
    ASSERT_DBG(global_enriched_dofs_[quant].size() > 0);
//     ASSERT_DBG(global_enriched_dofs_[quant][0].size() > 0);
    // if quantity not enriched, then returns 0
    unsigned int count;
    for(auto& q :global_enriched_dofs_)
        count += q.size();
    return count;
}

unsigned int XFEMElementDataBase::n_enriched_dofs(Quantity quant,
                                                  unsigned int local_enrichment_index) const
{
    ASSERT_DBG(quant < global_enriched_dofs_.size());
    ASSERT_DBG(local_enrichment_index < global_enriched_dofs_[quant].size());
    return global_enriched_dofs_[quant][local_enrichment_index].size();
}

unsigned int XFEMElementDataBase::n_enrichments_intersect() const
{
    unsigned int count = 0;
    for(unsigned int w=0; w < this->n_enrichments(); w++){
        if(enrichment_intersects_[w]) count++;
    }
    return count;
}


// void XFEMElementDataBase::print(ostream& out) const
// {
//     out << this << "xdata: ele " << this->ele_global_idx_;
//     out << " dofs[ ";
// //     for(unsigned int q=0; q<global_enriched_dofs_.size(); q++)
//     for(unsigned int q=0; q<2; q++)
//         for(unsigned int w=0; w<n_enrichments(); w++)
//             for(unsigned int j=0; j<global_enriched_dofs_[q][w].size(); j++){
//                 out << global_enriched_dofs_[q][w][j] << " ";
//             }
//     out << "]\n";
// }



template<int dim> XFEMElementSingularData<dim>::XFEMElementSingularData(){}
template<int dim> XFEMElementSingularData<dim>::~XFEMElementSingularData(){}


template<int dim>
void XFEMElementSingularData<dim>::create_sing_quads(ElementFullIter ele)
{
    ASSERT_DBG(dim == ele->dim());
//     ElementFullIter ele = mesh_->element(xdata.ele_global_idx());
    sing_quads_.resize(this->n_enrichments());
    this->enrichment_intersects_.resize(this->n_enrichments());
    
    MappingP1<dim,3> map;
    arma::mat proj = map.element_map(*ele);
    BoundingBox bb = ele->bounding_box();
    
    std::map<unsigned int, arma::vec> unit_points_inside;
    
    for(unsigned int w=0; w < this->n_enrichments(); w++){
        const Geometry& geom = this->enrichment_func(w)->geometry();
        unit_points_inside.clear();
        
//         DBGCOUT(<< "test q_points\n");
        for(unsigned int q=0; q < geom.q_points().size(); q++){
            const Space<3>::Point & p = geom.q_points()[q];
            
            //fast check
            if( ! bb.contains_point(p)) continue;
            
            //slow projection
            arma::vec unit_p = map.project_real_to_unit(p, proj);
            
            if( map.is_point_inside(unit_p)){
        
//                 DBGCOUT(<< "qpoint inside\n");
                unit_points_inside[q] = unit_p;
            }
        }
        
        QXFEM<dim,3>& qxfem = sing_quads_[w];
        qxfem.resize(unit_points_inside.size());
        double weight =  geom.effective_surface() / geom.q_points().size();
        std::map<unsigned int, arma::vec>::const_iterator pair;
        unsigned int i = 0;
        for(pair = unit_points_inside.begin(); pair != unit_points_inside.end(); pair++, i++){
            
            qxfem.set_point(i, RefElement<dim>::bary_to_local(pair->second));
            qxfem.set_real_point(i, geom.q_points()[pair->first]);
            qxfem.set_weight(i,weight);
        }
        //determine if enrichment is cross-secting the element
        this->enrichment_intersects_[w] = (qxfem.size() > 0);
        
        //debug output
        if(qxfem.size() > 0){
            DBGCOUT(<< "create_sing_quads on ele " << ele->index() << "\n");
            DBGCOUT(<< "quad[" << this->global_enrichment_index(w) << "] size " << sing_quads_[w].size() << "\n");
            DBGVAR(this->n_enrichments_intersect());
//             QXFEMFactory qfact;
//             qfact.gnuplot_refinement<dim>(ele,
//                                 FilePath("qxfem/", FilePath::output_file),
//                                 qxfem);
        }
    }
}


template<int dim>
const QXFEM<dim,3>& XFEMElementSingularData<dim>::sing_quadrature(unsigned int local_enrichment_index) const
{   ASSERT_DBG(local_enrichment_index < sing_quads_.size());
    return sing_quads_[local_enrichment_index];
}


// explicit instantiation
template class XFEMElementSingularData<2>;
template class XFEMElementSingularData<3>;
