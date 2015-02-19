/*
 * fe.hh
 *
 *  Created on: Dec 14, 2013
 *      Author: jb
 *
 *   TODO:
 *   - we need fast evaluation only for scalar valued Polynomial spaces
 *   - may be good to evaluate basis of several spaces together (keep monomial values,
 *     we can have some cache object for that
 *
 *   - it seems that we need not to evaluate vector or Tensor values (i.e. arma::vec<dim> and
 *     arma::mat_fixed<dim,dim> directly at least not for speed reasons
 *   - but it may be good for code quality
 *   - Think also about FEValues interface:
 *
 *     FEValuesSystem fe_values(quadrature, fe_system);
 *
 *     fe_values.reinit(dof_handler_cell)
 *
 *     fe_values[quantity].shape(i, q) // shape function i in quadrature point q
 *     fe_values[quantity] // has type  FEValues<ValueType>, that is just one quantity
 *                         // how to achieve this, quantitty has to have type ValueType
 *                         // but contain just index of the quantity in FeSystem
 *                         // to have proper resolution of [] operator for direrent
 *                         // return values
 *                         // precomputed values are saved in simple arrays
 *                         // ... be more precise in this
 *     value_extractor(fe_values, dof handler, vector) ... how to treat this
 *     ? how to treat FAValuesFaces
 *
 */

#ifndef FE_HH_
#define FE_HH_

/**
 * Base class for polynomial spaces, provides all common stuff safe
 * particular choice of monomial base.
 * TODO:
 * -
 */
template <int dim>
class PSpaceBase {
public:
	/**
	 * tensor product - put first space to the first coordinate and then others taking cartesian product..???
	 */
	PSpaceBase<dim> operator *(const PSpaceBase<dim> &other);
	/**
	 * direct addition - cartesian product
	 */
	PSpaceBase operator +(const PSpaceBase<dim> &other);
	/**
	 * Union - union of monomial sets
	 */
	PSpaceBase operator &(const PSpaceBase<dim> &other);


protected:
	/**
	 * mono_powers[i][j][k] is power of x_k coordinate in j-th component of  i-th "vector monomial",
	 * monomials are sorted
	 * according to ordering ... has to be determined  1,x,x^2, y, yx, yx^2, y^2, z, zx, zy, z^2
	 */
	vector<vector<vector<int>>> mono_powers_;

	unsigned int dimension_;

};

/********************
 * Elementary raw spaces
 * .. should be factory functions producing the space as instance of PSpaceBase
 */

/// One dimensional space of cooardinate "coord".
template <int dim>
class PSpaceX : public PSpaceBase<dim> {
public:
	PSpaceX(unsigned int coord) {
		dimension_=1;
		vector<vector<int>> one(1,vector<int>(dim));
		one[1][coord]=1;
		mono_powers_.push_back(one);
	}
};

///
template <int dim>


/**
 * Particular FE is given by N-dim monomial space and N basis polynomials given as
 * coefficients of N monomials.
 *
 * - what about vector valued FE
 * - how to compute values efficiently
 *   ... when used in XFEM the best approach is to compute all basis functions in one point at once
 *   then we compute values of monomials first and then their linear combinations as basis
 *
 *   Vector valued is more complex, since we have "vector valued" monomials as Raw space
 *   and then again their linear combinations are basis functions.
 */
class FiniteElement {

};



// Spaces form Brezzi Fortin
//
// H1, H2 conforming

// scalar polynomials with degree up to <p> degree
// BASIC space, need factory function
Pk<dim>(unsinged int degree)

// P kx,ky 2d - tensor product 
tensor_product(Pk<1>(kx), Pk<1>(ky))

// P kx,ky,kz 3d - tensor product 
tensor_product(Pk<1>(kx), Pk<1>(ky), Pk<1>(kz) )

// special case
Qk<1>(k) = Pk<1>(k)
Qk<2>(k) = tensor_product(Pk<1>(k), Pk<1>(k))
Qk<3>(k) = tensor_product(Pk<1>(k), Pk<1>(k))

// RT elements raw space
Pk^n + xPk
/*
... so we need:
- pk for all dimensions
- tensor product polynomial on R -> polynomial on R^n
- vector power scalar -> vector -> tensor valued
- vector X space (x,y,z) vector monomial
- addition - thats union of monomial basis
- product - is cartesian product of monomial sets

// Argyris triangle

spaces as subspaces with additional constrain
.. RT_k^0


*/
///////////////////////////////////////////////
// Boundary spaces

Rk ... Pk on sides

Tk ... Rk \cat continuous over boundary




#endif /* FE_HH_ */
