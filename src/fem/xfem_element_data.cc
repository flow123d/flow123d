#include "xfem_element_data.hh"
#include "singularity.hh"

#include <armadillo>

using namespace std;

XFEMElementSingularData::XFEMElementSingularData()
{}

XFEMElementSingularData::~XFEMElementSingularData()
{}


vector< vector< vector< int > > >& XFEMElementSingularData::global_enriched_dofs()
{
    return global_enriched_dofs_;
}


const vector< int >& XFEMElementSingularData::global_enriched_dofs(Quantity quant,
                                                                   unsigned int local_enrichment_index) const
{   ASSERT_DBG(quant < global_enriched_dofs_.size());
    ASSERT_DBG(local_enrichment_index < global_enriched_dofs_[quant].size());
    return global_enriched_dofs_[quant][local_enrichment_index];
}

unsigned int XFEMElementSingularData::n_enriched_dofs(Quantity quant) const
{
//     DBGVAR(ele_global_idx_); cout<<endl;
//     DBGVAR(global_enriched_dofs_.size()); cout<<endl;
//     DBGVAR(global_enriched_dofs_[quant].size()); cout<<endl;
//     DBGVAR(global_enriched_dofs_[quant][0].size()); cout<<endl;
    
    ASSERT_DBG(quant < global_enriched_dofs_.size());
    ASSERT_DBG(global_enriched_dofs_[quant].size() > 0);
//     ASSERT_DBG(global_enriched_dofs_[quant][0].size() > 0);
    // if quantity not enriched, then returns 0
    //FIXME: supposing that all enrichments have the same number of enr dofs
    return global_enriched_dofs_[quant].size() * global_enriched_dofs_[quant][0].size();
}

unsigned int XFEMElementSingularData::n_enriched_dofs(Quantity quant,
                                                      unsigned int local_enrichment_index) const
{
    ASSERT_DBG(quant < global_enriched_dofs_.size());
    ASSERT_DBG(local_enrichment_index < global_enriched_dofs_[quant].size());
    return global_enriched_dofs_[quant][local_enrichment_index].size();
}

void XFEMElementSingularData::create_sing_quads(ElementFullIter ele)
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
        double weight =  sing->circumference() / n_qpoints;
        std::map<unsigned int, arma::vec>::const_iterator pair;
        for(pair = unit_points_inside.begin(); pair != unit_points_inside.end(); pair++){
            
            qxfem.set_point(pair->first, RefElement<2>::bary_to_local(pair->second));
            qxfem.set_real_point(pair->first, sing->q_points()[pair->first]);
            qxfem.set_weight(pair->first,weight);
        }
        DBGCOUT(<< "quad[" << global_enrichment_index(w) << "] size " << sing_quads_[w].size() << "\n");
    }
}

unsigned int XFEMElementSingularData::n_singularities_inside() const
{
    unsigned int count = 0;
    for(unsigned int w=0; w < n_enrichments(); w++){
        if(is_singularity_inside(w)) count++;
    }
    return count;
}


bool XFEMElementSingularData::is_singularity_inside(unsigned int local_enrichment_index) const
{
    ASSERT_DBG(local_enrichment_index < sing_quads_.size());
    return sing_quads_[local_enrichment_index].size() > 0;
}


void XFEMElementSingularData::print(ostream& out) const
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


const QXFEM< int(2), int(3) >& XFEMElementSingularData::sing_quadrature(unsigned int local_enrichment_index) const
{   ASSERT_DBG(local_enrichment_index < sing_quads_.size());
    return sing_quads_[local_enrichment_index];
}
