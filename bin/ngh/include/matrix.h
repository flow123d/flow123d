#ifndef matrixH
#define matrixH

#include <iostream>
#include <armadillo>

typedef enum {
    no_solution,
    one_solution,
    inf_solutions,
    badconditioned,
    singular
} TNSolutions;

class TMatrix {
private:
    int nc;
    int nr;

    double* elm;

public:
    TMatrix(int);
    TMatrix(int, int);
    ~TMatrix();

    int NRows() const;
    int NCols() const;

    void Set(int, int, double);
    double Get(int, int) const;

    void SwapRows(int, int);

    friend std::ostream & operator <<(std::ostream&, const TMatrix&);
};

class TMVector {
private:
    int size;

    double* elm;

public:
    TMVector(int);
    ~TMVector();

    void Set(int, double);
    double Get(int);
    void SwapElements(int, int);

    friend std::ostream & operator <<(std::ostream&, const TMVector&);
};

TNSolutions Gauss(const TMatrix&, TMVector*, const TMVector&);

#endif
