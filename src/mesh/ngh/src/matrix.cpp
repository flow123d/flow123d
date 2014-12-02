#include <cmath>
#include "system/exc_common.hh"
#include "mesh/ngh/include/matrix.h"
#include "mesh/ngh/include/mathfce.h"

using namespace mathfce;

TMatrix::TMatrix(int size) {
    nc = size;
    nr = size;
    elm = new double[ nc * nr ];
}

TMatrix::TMatrix(int num_rows, int num_cols) {
    nc = num_cols;
    nr = num_rows;
    elm = new double[ nc * nr ];
}

TMatrix::TMatrix(const TMatrix & x)
{
    nc = x.nc;
    nr = x.nr;
    elm = new double[nc*nr];
    memcpy(elm,x.elm,nc*nr*sizeof(double));
}

TMatrix::~TMatrix() {
    delete[] elm;
}

TMVector::TMVector(int size) {
    this->size = size;
    elm = new double[size];
}

TMVector::TMVector(const TMVector & x)
{
    size=x.size;
    elm = new double[size];
    memcpy(elm,x.elm,size*sizeof(double));
}

TMVector::~TMVector() {
    delete[] elm;
}

std::ostream & operator <<(std::ostream& stream, const TMatrix& M) {
    for (int i = 1; i <= M.nr; i++) {
        for (int j = 1; j <= M.nc; j++) {
            stream << M.Get(i, j) << " ";
        }
        stream << "\n";
    }
    return stream;
}

std::ostream & operator <<(std::ostream& stream, const TMVector& V) {
    for (int i = 1; i <= V.size; i++) {
        stream << V.elm[ i - 1 ] << "\n";
    }
    return stream;
}

void TMatrix::Set(int row, int col, double value) {
    if (row > nr)
        THROW( ExcAssertMsg() << EI_Message( "Number of the row is greater than number of rows in the matrix.") );
    if (col > nc)
        THROW( ExcAssertMsg() << EI_Message( "Number of the column is greater than number of columns in the matrix.") );
    elm[ (row - 1) * nc + col - 1 ] = value;
    return;
}

double TMatrix::Get(int row, int col) const {
    if (row > nr)
        THROW( ExcAssertMsg() << EI_Message( "Number of the row is greater than number of rows in the matrix.") );
    if (col > nc)
        THROW( ExcAssertMsg() << EI_Message( "Number of the column is greater than number of columns in the matrix.") );
    return elm[ (row - 1) * nc + col - 1 ];
}

void TMatrix::SwapRows(int r1, int r2) {
    if (r1 > nr || r2 > nr) {
        THROW( ExcAssertMsg() << EI_Message( "Number of the row is greater than number of rows in the matrix.") );
    }

    for (int i = 1; i <= nc; i++) {
        double tmp = Get(r1, i);
        Set(r1, i, Get(r2, i));
        Set(r2, i, tmp);
    }

    return;
}

void TMVector::SwapElements(int i1, int i2) {
    if (i1 > size || i2 > size) {
        THROW( ExcAssertMsg() << EI_Message( "Number of the element is greater than size of the vector.") );
    }

    double tmp = elm[ i1 - 1 ];
    elm[ i1 - 1 ] = elm[ i2 - 1 ];
    elm[ i2 - 1 ] = tmp;

    return;
}

double TMVector::Get(int i) {
    if (i > size) {
    	THROW( ExcAssertMsg() << EI_Message("Number of the element is greater than size of the vector.") );
    }


    return elm[ i - 1 ];
}

void TMVector::Set(int i, double value) {
    if (i > size)
    	THROW( ExcAssertMsg() << EI_Message("Number of the element is greater than size of the vector.") );
    elm[ i - 1 ] = value;
    return;
}

int TMatrix::NRows() const {
    return nr;
}

int TMatrix::NCols() const {
    return nc;
}

TNSolutions Gauss(const TMatrix& A, TMVector* X, const TMVector& B) {
	TMatrix M(A);
    TMVector b(B);

    for (int i = 1; i < M.NRows(); i++) {
        double tmp = fabs(M.Get(i, i));
        int row = -1;
        for (int j = i + 1; j <= M.NRows(); j++)
            if (fabs(M.Get(j, i)) > tmp) {
                tmp = fabs(M.Get(j, i));
                row = j;
            }
        if (tmp < epsilon) {
            continue;
        }

        if (row != -1) {
            M.SwapRows(i, row);
            b.SwapElements(i, row);
        }

        for (int j = i + 1; j <= M.NRows(); j++) {
            tmp = M.Get(j, i) / M.Get(i, i);
            for (int k = i; k <= M.NCols(); k++) {
                M.Set(j, k, M.Get(j, k) - tmp * M.Get(i, k));
            }
            b.Set(j, b.Get(j) - tmp * b.Get(i));
        }
    }

    for (int i = M.NRows(); i >= 1; i--) {
        double tmp = b.Get(i);

        for (int k = i + 1; k <= M.NCols(); k++) {
            tmp -= M.Get(i, k) * X->Get(k);
        }

        if (IsZero(M.Get(i, i))) {
            if (IsZero(tmp)) {
                return inf_solutions;
            } else {
                return no_solution;
            }
        }

        X->Set(i, tmp / M.Get(i, i));
    }

    return one_solution;
}
