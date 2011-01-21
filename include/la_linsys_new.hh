/**
 * Classes for block linear systems. Namely it should provide following functionality:
 * 1) interface to various PETSC (or possibly other) matrices, with preallocation
 *    through the set_values method
 * 2) For coupling problems creation of the global linear system by connectiong sub systems of individual
 *    equations. Working cycle could be:
 *    1. creation of sub problems up to the global problem (just matrix dimensions and types)
 *    2. pre assembly - > preallocation
 *    3. assembly - by the same functions
 *    4. matrix modifications (schur complements, assembly of one global matrix, preconditioners ...
 *    5. iterative solution in parallel
 *    6. for nonlinear/time problems: jump to 3. and reuse the structures
 *
 * 3) Solution through explicit or implicit schur complements
 * 4) Preconditioning for coupled problems .
 * 5) In assembly, the access to the matrices in LinSys hierarchy is through successive call of block(i,j) method, but
 *    there also should be a class or another mechanism how to represent a coordinate in the hierarchy. Then one can access a MatriSimple
 *    directly by two coordinates and constantly access correct part of rhs vector
 *
 * Virtual test cases:
 * - Only water flow solved by 1, 2, 3, or 4 (Newton BC) Schur complements, explict and implicit
 *   Ideal: switch between Implicit and explicit only changing the particular class or some thing similarly simple
 * -
 *
 * Unsolved design problems:
 * - PETSC neumoznuje blokove vektory, tj. vektory by mely byt globalni a z nich by se mely pouzivat podvektory
 *   mozne reseni: blokovy globalni vektor se hierarchicky vytvari (pouze velikosti) a pak slouzi jako
 *   "fabrika" na globalni vektory i subvektory a nebo pro pristup k subvektorum globalniho vektoru
 *   .. detaily
 *
 */



/**
 * Abstract class for frontend for PETSc matrices.
 * purpose:
 * - allow preallocation through double assembly
 * - make interface to PETSC since PETSC evolve and it is simpler to make changes on one place
 * - make documented interface to PETSC and document used PETSC functionality this way
 *
 * Questions:
 *  - how to setup column distribution. It is necessary only for some of PETSC formats, isn't it?
 */
class MatrixSimple 
{public:
    /// possible states of the matrix
    typedef enum {
        INSERT=INSERT_VALUES,
        ADD=ADD_VALUES,
        ALLOCATE,
        DONE,
        NONE
    } SetValuesMode;

    /// Construct a parallel system with given local size.
    MatrixSimple(unsigned int loc_row_size, unsigned int col_size);

    /// multiply PETSC vector: y=Ax
    void multiply(Vec x, Vec y);
    /// multiply vector and ad onother one: y=Ax+b
    void multiply_add(Vec x, Vec b, Vec y);
    /// multiply matrix .. has to deal with preallocation of result and reuse of structure
    void multiply_mat(...);
    /// add and scale sparse matrices ... has to deal with different nonzero structures,
    /// in such a case petsc rutines are inceribly slow
    void axpy(...);
    void aypx(...);

    /// @name access members @{
    /// Access to distribution of rows, thhrou that you can get also matix size and local size. But for convenience
    /// we provide also direct functions
    Distribution &get_row_ds();
    Distribution &get_col_ds();
    /// Get global system size.
    inline unsigned int size();
    /// Get matrix. SHOULD NOT BE USED !!!
    inline const Mat &get_matrix()
    { return matrix; }
    /// @}

    virtual void start_allocation()=0;
    void start_add_assembly();
    void start_insert_assembly();
    void finalize(MatAssemblyType assembly_type=MAT_FINAL_ASSEMBLY);

    virtual void preallocate_matrix()=0;
    virtual void preallocate_values(int nrow,int *rows,int ncol,int *cols)=0;

    virtual void view_local_matrix()=0;

    /// Set full rectangular submatrix of the system matrix.
    void mat_set_values(int nrow,int *rows,int ncol,int *cols,PetscScalar *vals);

    /// Set one element of the system matrix.
    inline void mat_set_value(int row,int col,PetscScalar val);

    inline void set_symmetric(bool flag = true);
    inline bool is_symmetric();

    inline void set_positive_definite(bool flag = true);
    inline bool is_positive_definite();

    /// Output the matrix in the Matlab format possibly with given renumbering of rows/cols
    void view(std::ostream output_stream, int * row_mapping = NULL, int * col_mapping = NULL);

    virtual ~LinSys();

protected:
    Distribution row_ds;            ///< Distribution of continuous blocks of system rows among the processors.
    Distribution col_ds;            ///< Distribution of multiplied vector.
    bool     symmetric;             ///< Flag for the symmetric system.
    bool     positive_definite;     ///< Flag for positive definite system.
    SetValuesMode status;        ///< Set value status of the linear system.

    Mat     matrix;                  ///< Petsc matrix of the problem.

};

/**
 *  Array of MatrixSimple
 *
 */
class MatrixArray 
{
private:
  std::vector<std::vector<MatrixSimple *>> mat_array;

public:
  MatrixArray ... vice konstruktoru pro ruzne slepovani MatrixSimple a jiz existujicich MatrixArray
  
  MatrixSimple *block(i_row, i_col)
  multiply_vec( Vector x, Vector y)
}

/**
 * Array or hierarchy over the global vector.
 */
class VecArray
{

};


/**
 * LinSys - matrix and a particular way how to compute solution for given RHS, i.e. this class can perform
 * action of  the matrix inverse
 *
 * This is abstract class for members of possible hierarchical tree of the whole system.
 *
 */
class LinSys {
  LinSys *block(i,j)
  MatrixArray *matrix()

  virtual solve( Vector solution, Vector RHS);
}

/**
 * The LinSys for direct factorization, its matrix has to be simple. Factorization will be made through PETSc
 * so the used backend library library can be influenced by on-line parameters.
 */
class LinSysDirect : public LinSys 
{
    MatrixSimple *block()
}

/**
 * Array of four submatrices A00 A10 A01 A11 where A00 has to be LinSys and A11 has to be Simple.
 * Solves the system by recursive solution of schur complement system. Can be implemented as
 * interface to PETSC functionality see PCFIELDSPLIT
 */
class LinSysImplicitSchur : public LinSys 
{
    LinSys *block(i,j)
    MatrixArray *matrix()
}
/**
 * This LinSys should after assembly construct an explicit Schur complemnt system. User has to provide
 * an LinSysDirect for one of the diagonal blocks.
 *
 * There is a problem how to make possible several successive schur complements.
 * This can not be organised bottom up like ImplicitSchur.
 */
class LinSysExplicitSchur : public LinSys
{
}

/**
 * Array of four matrices 00 01 10 11, where 00 and 11 are LinSys and their application is used as preconditioner for the whole system.
 * again see PCFIELDSPLIT
 */
class LinSysJacobiCoupling : public LinSys 
{
    LinSys *block11()
    MatrixArray *block12()
    MatrixArray *block21()
    LinSys *block22()
}





