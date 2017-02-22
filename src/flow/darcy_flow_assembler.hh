


#ifndef DARCY_FLOW_ASSEMBLER_HH_
#define DARCY_FLOW_ASSEMBLER_HH_

#include "darcy_flow_mh.hh"
#include "richards_lmh.hh"
#include "darcy_flow_assembly.hh"
#include "assembly_lmh.hh"


class AssemblerBase
{
public:
    typedef std::shared_ptr<DarcyMH::EqData> AssemblyDataPtr;
    typedef std::vector<std::shared_ptr<AssemblyBase> > MultidimAssembly;

    AssemblerBase(AssemblyDataPtr d)
    : ad_(d), local_boundary_index(0)
    {
        for(int i=0; i<1000; i++) zeros[i]=0.0;
    }
    
    virtual ~AssemblerBase(){}
    
    AssemblyDataPtr ad_;
    MultidimAssembly dim_assembler;
    
    void assemble(LocalElementAccessorBase<3> ele_ac){
        
        unsigned int dim = ele_ac.dim();
        dim_assembler[dim-1]->assemble(ele_ac);
        
        //assembly_dim_connections(ele_ac);
        
        schur_allocations(ele_ac);
        
        if (ad_->balance != nullptr)
            add_fluxes_in_balance_matrix(ele_ac);
    }
    
protected:
    /// Auxiliary counter for boundary balance.
    unsigned int local_boundary_index;
    
    int tmp_rows[100];
    // to make space for second schur complement, max. 10 neighbour edges of one el.
    double zeros[1000];
    
    int rows[2];
    double local_vb[4]; // 2x2 matrix
    int edge_rows[4];
    
    void add_fluxes_in_balance_matrix(LocalElementAccessorBase<3> ele_ac){
        
        for (unsigned int i = 0; i < ele_ac.n_sides(); i++) {
            Boundary* bcd = ele_ac.side(i)->cond();

            if (bcd) {
                /*
                    DebugOut()("add_flux: {} {} {} {}\n",
                            mh_dh.el_ds->myp(),
                            ele_ac.ele_global_idx(),
                            local_boundary_index,
                            side_row);*/
                ad_->balance->add_flux_matrix_values(ad_->water_balance_idx, local_boundary_index,
                                                     {ele_ac.side_row(i)}, {1});
                ++local_boundary_index;
            }
        }
    }

    void schur_allocations(LocalElementAccessorBase<3> ele_ac){
            
        LinSys* ls = ad_->lin_sys;
        Neighbour *ngh;
        int ele_row = ele_ac.ele_row();
        unsigned int nsides = ele_ac.n_sides();
        
        // edge dofs are saved in local system
//             LocalSystem& loc = dim_assembler[ele_ac.dim()-1]->get_local_system();
//             int* edge_rows = & loc.row_dofs[ele_ac.dim() + 2];
        for (unsigned int i = 0; i < nsides; i++)
            edge_rows[i] = ele_ac.edge_row(i);
        
        for (unsigned int i = 0; i < ele_ac.full_iter()->n_neighs_vb; i++) {
            // every compatible connection adds a 2x2 matrix involving
            // current element pressure  and a connected edge pressure
            ngh= ele_ac.full_iter()->neigh_vb[i];
            tmp_rows[0]=ele_row;
            tmp_rows[1]=ad_->mh_dh->row_4_edge[ ngh->edge_idx() ];

            if (ad_->n_schur_compls == 2) {
                // for 2. Schur: N dim edge is conected with N dim element =>
                // there are nz between N dim edge and N-1 dim edges of the element

                ls->mat_set_values(nsides, edge_rows, 1, tmp_rows+1, zeros);
                ls->mat_set_values(1, tmp_rows+1, nsides, edge_rows, zeros);

                // save all global edge indices to higher positions
                tmp_rows[2+i] = tmp_rows[1];
            }
        }

        // add virtual values for schur complement allocation
        uint n_neigh;
        switch (ad_->n_schur_compls) {
        case 2:
            n_neigh = ele_ac.full_iter()->n_neighs_vb;
            // Connections between edges of N+1 dim. elements neighboring with actual N dim element 'ele'
            OLD_ASSERT(n_neigh*n_neigh<1000, "Too many values in E block.");
            ls->mat_set_values(ele_ac.full_iter()->n_neighs_vb, tmp_rows+2,
                    ele_ac.full_iter()->n_neighs_vb, tmp_rows+2, zeros);

        case 1: // included also for case 2
            // -(C')*(A-)*B block and its transpose conect edge with its elements
            ls->mat_set_values(1, &ele_row, nsides, edge_rows, zeros);
            ls->mat_set_values(nsides, edge_rows, 1, &ele_row, zeros);
            // -(C')*(A-)*C block conect all edges of every element
            ls->mat_set_values(nsides, edge_rows, nsides, edge_rows, zeros);
        }
    }

};

class AssemblerMH : public AssemblerBase
{
public:
    AssemblerMH(AssemblyDataPtr d)
    : AssemblerBase(d)
    {
        dim_assembler = AssemblyBase::create< AssemblyMH >(d);
    }
};

class AssemblerLMH : public AssemblerBase
{
public:
    AssemblerLMH(std::shared_ptr<RichardsLMH::EqData> d)
    : AssemblerBase(d)
    {
        dim_assembler = AssemblyBase::create< AssemblyLMH >(d);
    }
};


#endif // DARCY_FLOW_ASSEMBLER_HH_
