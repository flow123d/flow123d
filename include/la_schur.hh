/**
 *
 * $Id: big_la_schur_complement.hh$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief  Assembly explicit Schur complement for the given linear system.
 * Provides method for resolution of the full original vector of unknowns.
 *
 */

#ifndef LA_SCHUR_HH_
#define LA_SCHUR_HH_

struct Solver;
class LinSys;

/**
 * @brief Schur complement class for a PETSC based linear system
 */

typedef enum SchurState {
    created,    // created or after reuse
    formed,     // formed complement system redy to solve
    solved      // solved ready to resolve
} SchurState;

typedef class SchurComplement {
public:
    SchurComplement(LinSys *orig,Mat inv_a, IS ia = NULL);

    LinSys *get_system() const {return (Compl);}
    LinSys *get_orig_system() const {return (Orig);}
    Mat get_a_inv() const {return (IA);}
    void set_spd();
    void reuse() {state=created;}

    void scale(double factor);
    void solve(Solver *solver);
    ~SchurComplement();

    // TODO: should be at least protected
    void form_schur();
    void resolve();

private:
    Mat IA;                     // Inverse of block A
    Mat B, Bt;                   // B and B' block (could be different from real B transpose)
    Mat xA;                     // Bt*IA*B
    Mat IAB;                    // reconstruction matrix IA * B
    int locSizeA, locSizeB;     // loc size of the A and B block
    IS IsA, IsB;                // parallel index sets of the A and B block
    IS fullIsA,fullIsB;         // whole IsA  and IsB on each proc
    Vec RHS1, RHS2;             // A and B - part of the RHS
    Vec Sol1, Sol2;             // A and B part of solution
    MatReuse mat_reuse;        // reuse structures after first computation of schur
    SchurState state;           // object internal state

                                //                A  B     Sol1      RHS1
    LinSys *Orig;     // Original Linear System:  B' C  *  Sol2  =   RHS2
    LinSys *Compl;    // Schur complement system: (C - B' IA B) * Sol2 = (B' * IA * RHS1 - RHS2)
} SchurComplement;

#endif /* LA_SCHUR_HH_ */
