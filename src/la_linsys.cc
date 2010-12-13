/*!
 * $Id: la_linsys.cc 613 2010-07-27 07:26:11Z jan.brezina $
 * $Revision: 613 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2010-07-27 09:26:11 +0200 (Tue, 27 Jul 2010) $
 *
 * @brief   Wrappers for linear systems based on MPIAIJ and MATIS format.
 * @author  Jan Brezina
 *
 */

#include <algorithm>
#include <petscmat.h>
#include <system.hh>
#include <la_linsys.hh>

/**
 *  @brief Constructs a parallel system with given local size.
 *
 *  By the parameter @param vec_lsize we initialize distribution of the rhs and solution vectors.
 *  For MPIAIJ matrix this is also distribution of its rows, but for IS matrix this is only size of
 *  principial local part without interface.
 */

LinSys::LinSys(unsigned int vec_lsize, double *sol_array)
:vec_ds(vec_lsize),symmetric(false),positive_definite(false),status(NONE)
{
    // create PETSC vectors
    v_rhs=(double *) xmalloc(sizeof(double) * (this->vec_lsize() + 1) );
    VecCreateMPIWithArray(PETSC_COMM_WORLD, this->vec_lsize(), PETSC_DECIDE, v_rhs, &rhs);
    VecZeroEntries(rhs);

    if (sol_array == NULL) v_solution=(double *)xmalloc(sizeof(double) * (this->vec_lsize() + 1));
    else v_solution=sol_array;
    VecCreateMPIWithArray(PETSC_COMM_WORLD, this->vec_lsize(), PETSC_DECIDE, v_solution, &solution);
    own_solution=false;
    //VecZeroEntries(solution);
}


void LinSys::start_add_assembly()
{
    switch (status) {
        case ALLOCATE:
            preallocate_matrix();
            break;
        case INSERT:
            finalize(MAT_FLUSH_ASSEMBLY);
            break;
        case ADD:
        case DONE:
            break;
        default:
            ASSERT(0, "Can not set values. Matrix is not preallocated.\n");
    }
    status=ADD;
}

void LinSys::start_insert_assembly()
{
    switch (status) {
        case ALLOCATE:
            preallocate_matrix();
            break;
        case ADD:
            finalize(MAT_FLUSH_ASSEMBLY);
            break;
        case INSERT:
        case DONE:
            break;
        default:
            ASSERT(0, "Can not set values. Matrix is not preallocated.\n");
    }
    status=INSERT;
}


void LinSys::finalize(MatAssemblyType assembly_type)
{
    if (status == ALLOCATE) {
        xprintf(Warn, "Finalizing linear system without setting values.\n");
        preallocate_matrix();
    }
    MatAssemblyBegin(matrix, assembly_type);
    VecAssemblyBegin(rhs);
    if (assembly_type == MAT_FINAL_ASSEMBLY) status=DONE;
    MatAssemblyEnd(matrix, assembly_type);
    VecAssemblyEnd(rhs);
}

void view(std::ostream output_stream, int * output_mapping = NULL)
{

}

LinSys:: ~LinSys()
{
    MatDestroy(matrix);
    VecDestroy(rhs);
    VecDestroy(solution);

    xfree(v_rhs);
    if (own_solution) xfree(v_solution);
}

#if 0

// ======================================================================================
// LSView - output assembled system in side,el,edge ordering

void LSView(LinSystem *ls)
{
    FILE *f;
    int i, j, n;
    PetscInt ncols;
    const PetscInt *cols;
    const PetscScalar *vals;
    MatInfo info;

    F_ENTRY;

    // output matrix
    MatGetInfo(ls->A,MAT_GLOBAL_SUM,&info);
    MPI_Barrier(PETSC_COMM_WORLD);
    if (ls->ds->myp() == 0) {
        f=xfopen("matrix.dat","w");
        xfprintf(f,"%d %d 0.0\n",ls->size, ls->size);
        //xfprintf(f,"zzz=zeros(%d,3)\n",(int)(info.nz_used));
        //xfprintf(f,"zzz=[\n");
        xfclose(f);
    }
    for(n=0;n<ls->ds->np();n++) {
        MPI_Barrier(PETSC_COMM_WORLD);
        if (n==ls->ds->myp()) {
            f=xfopen("matrix.dat","a");
            for(i=ls->ds->begin(); i<ls->ds->end();i++) {
                MatGetRow(ls->A,i,&ncols,&cols,&vals);
                for(j=0;j<ncols;j++)
                    xfprintf(f,"%d %d %f\n",ls->old_4_new[i]+1,ls->old_4_new[cols[j]]+1,vals[j]);
                MatRestoreRow(ls->A,i,&ncols,&cols,&vals);
            }
            xfclose(f);
        }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    if (ls->ds->myp() == 0) {
        f=xfopen("rhs.m","w");
        xfprintf(f,"yyy = zeros(%d,2)\n yyy=[\n",ls->size);
        xfclose(f);
    }
    // output vec
    for(n=0; n<ls->ds->np(); n++) {
        MPI_Barrier(PETSC_COMM_WORLD);
        if (n==ls->ds->myp()) {
            f=xfopen("rhs.m","a");
            for(i=ls->ds->begin(); i<ls->ds->end(); i++) {
                xfprintf(f,"%d %f\n",ls->old_4_new[i]+1,ls->vb[i-ls->ds->begin()]);
            }
            xfclose(f);
        }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    if (ls->ds->myp() == 0) {
        f=xfopen("rhs.m","a");
        xfprintf(f,"];\n xxx=sortrows(yyy); Vec=xxx(:,2);\n");
        xfclose(f);
    }
}

// ======================================================================================
// LSView - output assembled system in side,el,edge ordering

void MyVecView(Vec v,int *order,const char* fname) {

FILE    *f;


int i,n;
double *array;


    F_ENTRY;
    Distribution ds(v); // get distribution of the vector


    if (ds.myp() == 0) {
        f=xfopen(fname,"w");
        xfprintf(f,"yyy = zeros(%d,2)\n yyy=[\n",ds.size());
        xfclose(f);
    }
    // output vec
    for(n=0; n<ds.np(); n++) {
        MPI_Barrier(PETSC_COMM_WORLD);
        if (n==ds.myp()) {
            VecGetArray(v,&array);
            f=xfopen(fname,"a");
            for(i=ds.begin(); i<ds.end(); i++) {
                xfprintf(f,"%d %f\n",order[i]+1,array[i-ds.begin()]);
            }
            VecRestoreArray(v,&array);
            xfclose(f);
        }
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    if (ds.myp() == 0) {
        f=xfopen(fname,"a");
        xfprintf(f,"];\n xxx=sortrows(yyy); Vec=xxx(:,2);\n");
        xfclose(f);
    }
}


//=========================================================================================
/*! @brief convert linear system to pure CSR format
 *
 *  - allocate CSR arrays according to MatGetInfo
 *  - indexes are form 0
 *  - make CSR format of the matrix (i.e. I,J,A vectors) from the PETSC format
 */

void LSSetCSR( LinSystem *mtx )
{
    int i,pos,row,nnz;
    PetscInt nnz_loc;
    const PetscInt *cols;
    const PetscScalar *vals;
    MatInfo info;

    F_ENTRY;

    MatGetInfo(mtx->A,MAT_LOCAL,&info);
    nnz=(int)(info.nz_used)+1; // allocate one more to be sure
    // allocate
    if (mtx->i == NULL) mtx->i=(int *)xmalloc( (mtx->size+1)*sizeof(int) );
    if (mtx->j == NULL) mtx->j=(int *)xmalloc( nnz*sizeof(int) );
    if (mtx->a == NULL) mtx->a=(double *)xmalloc( nnz*sizeof(double) );
    DBGMSG("mtx size: %d nnz: %d\n",mtx->size,nnz);
    // setup matrix from PETSC SeqAIJ format
    pos=0;
    for( row=0; row< mtx->size; row++ ) {
        mtx->i[row]=pos;
        MatGetRow(mtx->A,row,&nnz_loc,&cols,&vals);
        DBGMSG("CSR row: %d nnz_loc %d\n",row,nnz_loc);
        for(i=0; i<nnz_loc; i++) {
                ASSERT(pos<nnz,"More nonzeroes then allocated! row: %d entry: %d\n",row,i);
                mtx->j[pos]=cols[i];
                mtx->a[pos]=vals[i];
                pos++;
        }
        MatRestoreRow(mtx->A,row,&nnz_loc,&cols,&vals);
    }
    mtx->i[row]=pos;
}

void LSFreeCSR( LinSystem *mtx )
{
    xfree(mtx->i);
    xfree(mtx->j);
    xfree(mtx->a);
}

#endif

//**********************************************************************************************


void LinSys_MPIAIJ::start_allocation()
{
     if (status != NONE) {
         // reinit linear system

     }
     VecCreateMPI(PETSC_COMM_WORLD, vec_ds.lsize(), PETSC_DECIDE, &(on_vec));
     VecDuplicate(on_vec,&(off_vec));
     status=ALLOCATE;

     DBGMSG("allocation started\n");
}

void LinSys_MPIAIJ::preallocate_matrix()
{
     ASSERT(status == ALLOCATE, "Linear system has to be in ALLOCATE status.");

     PetscScalar *on_array, *off_array;
     int *on_nz, *off_nz;
     int i;

     // assembly and get values from counting vectors, destroy them
     VecAssemblyBegin(on_vec);
     VecAssemblyBegin(off_vec);

     on_nz=(int *)xmalloc( 2 * vec_ds.lsize() * sizeof(int));
     off_nz=on_nz + vec_ds.lsize();

     VecAssemblyEnd(on_vec);
     VecAssemblyEnd(off_vec);

     VecGetArray(on_vec,&on_array);
     VecGetArray(off_vec,&off_array);

     for(i=0; i<vec_ds.lsize(); i++) {
         on_nz[i]=(int)(on_array[i]+0.1);        // small fraction to ensure correct rounding
         off_nz[i]=(int)(off_array[i]+0.1);
     }

     VecRestoreArray(on_vec,&on_array);
     VecRestoreArray(off_vec,&off_array);
     VecDestroy(on_vec);
     VecDestroy(off_vec);

     // create PETSC matrix with preallocation
     MatCreateMPIAIJ(PETSC_COMM_WORLD, vec_ds.lsize(), vec_ds.lsize(), PETSC_DETERMINE, PETSC_DETERMINE,
              PETSC_NULL, on_nz, PETSC_NULL, off_nz, &matrix);

     if (symmetric) MatSetOption(matrix, MAT_SYMMETRIC, PETSC_TRUE);
     MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);


     xfree(on_nz);
}

void LinSys_MPIAIJ::preallocate_values(int nrow,int *rows,int ncol,int *cols)
{
     int i,j,row,col;

     for (i=0; i<nrow; i++) {
         row=rows[i];
         for(j=0; j<ncol; j++) {
             col=cols[j];
             if (vec_ds.get_proc(row) == vec_ds.get_proc(col))
                 VecSetValue(on_vec,row,1.0,ADD_VALUES);
             else
                 VecSetValue(off_vec,row,1.0,ADD_VALUES);
         }
     }
}

void LinSys_MPIAIJ::view_local_matrix()
{
     // print local subdomain matrix
     xprintf(Msg,"Printing of local matrix is not supported yet for MPIAIJ matrix. \n");

}


LinSys_MPIAIJ:: ~LinSys_MPIAIJ()
{
     if (status == ALLOCATE) {
         VecDestroy(on_vec);
         VecDestroy(off_vec);
     }
}

//**********************************************************************************************

LinSys_MATIS::LinSys_MATIS(unsigned int vec_lsize,  int subdomain_size, int *global_row_4_sub_row, double *sol_array)
: LinSys(vec_lsize, sol_array), subdomain_size(subdomain_size)
{
    PetscErrorCode err;

    // vytvorit mapping v PETSc z global_row_4_sub_row
    err = ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, subdomain_size, global_row_4_sub_row, &map_local_to_global);
    ASSERT(err == 0,"Error in ISLocalToGlobalMappingCreate.");
    xprintf(Msg,"Error code in ISLocalToGlobalMappingCreate %d \n",err);

    // initialize loc_rows array
    loc_rows_size=100;
    loc_rows = new int[loc_rows_size];

};


void LinSys_MATIS::start_allocation()
{
     PetscErrorCode err;

     if (status != NONE) {
         // reinit linear system

     }
     err = MatCreateIS(PETSC_COMM_WORLD,  vec_ds.lsize(), vec_ds.lsize(), vec_ds.size(), vec_ds.size(), map_local_to_global, &matrix);
     ASSERT(err == 0,"Error in MatCreateIS.");
     xprintf(Msg,"Error code MatCreateIS %d \n",err);

     err = MatISGetLocalMat(matrix, &local_matrix);
     ASSERT(err == 0,"Error in MatISGetLocalMat.");
     xprintf(Msg,"Error code MatISGetLocalMat %d \n",err);

     subdomain_nz= new int[subdomain_size];      // count local nozero for every row of subdomain matrix
     SET_ARRAY_ZERO(subdomain_nz,subdomain_size); // set zeros to the array

     status=ALLOCATE;

     DBGMSG("allocation started\n");
}

void LinSys_MATIS::preallocate_matrix()
{
     ASSERT(status == ALLOCATE, "Linear system has to be in ALLOCATE status.");

     // printing subdomain_nz
     DBGPRINT_INT("subdomain_nz",subdomain_size,subdomain_nz);
      

     // preallocation of local subdomain matrix
     MatSeqAIJSetPreallocation(local_matrix, PETSC_NULL, subdomain_nz);

     if (symmetric) MatSetOption(matrix, MAT_SYMMETRIC, PETSC_TRUE);
     MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

     delete[] subdomain_nz;
}

void LinSys_MATIS::preallocate_values(int nrow,int *rows,int ncol,int *cols)
{
     int i,row, n_loc_rows;
     int irow, indrow, indrow_loc;
     PetscErrorCode err;

     if (loc_rows_size < nrow) {
         delete[] loc_rows;
         loc_rows_size=nrow;
         loc_rows=new int[loc_rows_size];

     }

     // translate from global to local indexes
     err = ISGlobalToLocalMappingApply(map_local_to_global, IS_GTOLM_DROP, nrow, rows, &n_loc_rows, loc_rows);
     ASSERT(err == 0,"Error in ISGlobalToLocalMappingApply.");
     ASSERT(nrow == n_loc_rows,"Not all global indices translated to local indices.");
     //xprintf(Msg,"Error code ISGlobalToLocalMappingApply %d \n",err);
     // printing subdomain_embedding
     //DBGPRINT_INT("embed_element_to",nrow,loc_rows);


     /* We don't need columns for preallocation.

     if (loc_cols_size < ncol) {
              delete loc_cols;
              loc_cols_size=ncol;
              loc_cols=new int[loc_cols_size];
          }

     // translate from global to local indexes
     ISGlobalToLocalMappingApply(local_to_global, IS_GTOLM_DROP, ncol, cols, &n_loc_cols, loc_cols);
     */
     for (i=0; i<n_loc_rows; i++) {
         row=loc_rows[i];
         subdomain_nz[row] = subdomain_nz[row] + ncol;
     }
}

void LinSys_MATIS::view_local_matrix()
{
     PetscErrorCode err;
     PetscViewer lab;
     char numstring[6] = "00000";
     int myid;
     int ierr;

//     err = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF,PETSC_VIEWER_ASCII_DENSE);
//     ASSERT(err == 0,"Error in PetscViewerSetFormat.");

     // print local subdomain matrix
     err = MatView(local_matrix,PETSC_VIEWER_STDOUT_SELF);
     ASSERT(err == 0,"Error in MatView.");
//
//     ierr=MPI_Comm_rank(PETSC_COMM_WORLD, &(myid));
//     ASSERT( ! ierr , "Can not get MPI rank.\n" );
//
//     sprintf(numstring,"%5.5d",myid);
//     printf("Nunstring is >>>%s<<<\n",numstring);
//
//     PetscViewerASCIIOpen(PETSC_COMM_SELF,numstring,&lab);
//
//     ASSERT(!(loc_rows == NULL),"Local matrix is not assigned.");
//     err = PetscViewerSetFormat(lab,PETSC_VIEWER_ASCII_DENSE);
//     ASSERT(err == 0,"Error in PetscViewerSetFormat.");
//
//     // print local subdomain matrix
//     err = MatView(local_matrix,lab);
//     ASSERT(err == 0,"Error in MatView.");
//
//     PetscViewerDestroy(lab);
}

LinSys_MATIS:: ~LinSys_MATIS()
{
     PetscErrorCode err;

     // destroy local subdomain matrix
     err = MatDestroy(local_matrix);
     ASSERT(err == 0,"Error in MatDestroy.");
     xprintf(Msg,"Error code MatDestroy %d \n",err);

     // destroy mapping
     err = ISLocalToGlobalMappingDestroy(map_local_to_global);
     ASSERT(err == 0,"Error in ISLocalToGlobalMappingDestroy.");
     xprintf(Msg,"Error code %d \n",err);


     delete[] loc_rows;
     if (status == ALLOCATE) {
         delete subdomain_nz;
     }

}
