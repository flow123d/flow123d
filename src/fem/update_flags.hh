/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    update_flags.hh
 * @brief   Enum type UpdateFlags indicates which quantities are to be
 *          recomputed on each finite element cell.
 * @author  Jan Stebel
 */

#ifndef UPDATE_FLAGS_HH_
#define UPDATE_FLAGS_HH_

/**
 * \enum UpdateFlags
 * @brief Enum type UpdateFlags indicates which quantities are to be
 * recomputed on each finite element cell.
 *
 * Selecting these flags in a restrictive way is crucial for the
 * efficiency of FEValues::reinit() and FESideValues::reinit().
 * Therefore, only the flags actually
 * needed should be selected. It is the responsibility of the involved
 * Mapping and FiniteElement to add additional flags according to
 * their own requirements.

 * By default, all flags are off, i.e. no reinitialization will be
 * done.
 *
 * You can select more than one flag by concatenation
 * using the bitwise or operator|(UpdateFlags,UpdateFlags).
 *
 * <h3>Generating the actual flags</h3>
 *
 * When given a set of UpdateFlags @p flags, the FEValues object must
 * determine, which values will have to be computed once only for the
 * reference cell and which values will have to be updated for each
 * cell. Here, it is important to note that in many cases, the
 * FiniteElement will require additional updates from the Mapping. To
 * this end, the following auxiliary functions have been implemented:
 *
  * FiniteElement::update_each() determine the values required by
 * the FiniteElement on each cell. The same function exists in Mapping.
 *
 * FEValuesBase::update_each() is used to compute the union
 * of all values to be computed ever. It does this by first adding to
 * the flags set by the user all flags added by the FiniteElement.
 * This new set of flags is then given to the Mapping and all flags
 * required there are added.
 *
 * This union of all flags is given to Mapping::fill_fe_values() and
 * FiniteElement::fill_fe_values(), where the quantities indicated by
 * the flags are computed.
 *
 * The flags finally stored in FEValues then are the union of all the
 * flags required by the user, by FiniteElement and by Mapping, for
 * computation once or on each cell.
 */
enum UpdateFlags
{
                       //! No update
      update_default                      = 0,
                       //! Shape function values
                       /**
                    * Compute the values of the
                    * shape functions at the
                    * quadrature points on the
                    * real space cell. For the
                    * usual Lagrange elements,
                    * these values are equal to
                    * the values of the shape
                    * functions at the quadrature
                    * points on the unit cell, but
                    * they are different for more
                    * complicated elements, such
                    * as Raviart-Thomas
                    * elements.
                    */
      update_values                       = 0x0001,
                       //! Shape function gradients
                       /**
                    * Compute the gradients of the
                    * shape functions in
                    * coordinates of the real
                    * cell.
                    */
      update_gradients                    = 0x0002,
                       //! Transformed quadrature points
                       /**
                    * Compute the quadrature
                    * points transformed into real
                    * cell coordinates.
                    */
      update_quadrature_points            = 0x0004,
                       //! Transformed quadrature weights
                       /**
                    * Compute the quadrature
                    * weights on the real cell,
                    * i.e. the weights of the
                    * quadrature rule multiplied
                    * with the determinant of the
                    * Jacobian of the
                    * transformation from
                    * reference to real cell.
                    */
      update_JxW_values                   = 0x0008,
                                       //! Normal vectors
                       /**
                    * Compute the normal vectors,
                    * either for a face or for a
                    * cell of codimension
                    * one. Setting this flag for
                    * any other object will raise
                    * an error.
                    */
      update_normal_vectors               = 0x0010,
                       //! Volume element
                       /**
                    * Compute the Jacobian of the
                    * transformation from the
                    * reference cell to the real
                    * cell.
                    */
      update_jacobians                    = 0x0020,
                       //! Volume element
                       /**
                    * Compute the inverse
                        * Jacobian of the
                    * transformation from the
                    * reference cell to the real
                    * cell.
                    */
      update_inverse_jacobians            = 0x0040,
                       //! Determinant of the Jacobian
                       /**
                    * Compute the volume element
                    * in each quadrature point.
                    */
      update_volume_elements              = 0x0080,
                       //! Transformed quadrature weight for cell sides
                       /**
                    * Same as update_JxW_values but for quadratures living
                    * on a side of the cell.
                    */
      update_side_JxW_values              = 0x0100
};


/**
 * Output operator which outputs update flags as a set of or'd text values.
 *
 * @ref UpdateFlags
 */
template <class STREAM>
inline
STREAM& operator << (STREAM& s, UpdateFlags u)
{
  s << " UpdateFlags|";
  if (u & update_values)                       s << "values|";
  if (u & update_gradients)                    s << "gradients|";
  if (u & update_quadrature_points)            s << "quadrature_points|";
  if (u & update_JxW_values)                   s << "JxW_values|";
  if (u & update_normal_vectors)               s << "normal_vectors|";
  if (u & update_jacobians)                    s << "jacobians|";
  if (u & update_inverse_jacobians)            s << "inverse_jacobians|";
  return s;
}


/**
 * Global operator which returns an object in which all bits are set
 * which are either set in the first or the second argument. This
 * operator exists since if it did not then the result of the bit-or
 * <tt>operator |</tt> would be an integer which would in turn trigger
 * a compiler warning when we tried to assign it to an object of type
 * UpdateFlags.
 *
 * @ref UpdateFlags
 */
inline
UpdateFlags
operator | (UpdateFlags f1, UpdateFlags f2)
{
  return static_cast<UpdateFlags> (
    static_cast<unsigned int> (f1) |
    static_cast<unsigned int> (f2));
}




/**
 * Global operator which sets the bits from the second argument also
 * in the first one.
 *
 * @ref UpdateFlags
 */
inline
UpdateFlags &
operator |= (UpdateFlags &f1, UpdateFlags f2)
{
  f1 = f1 | f2;
  return f1;
}


/**
 * Global operator which returns an object in which all bits are set
 * which are set in the first as well as the second argument. This
 * operator exists since if it did not then the result of the bit-and
 * <tt>operator &</tt> would be an integer which would in turn trigger
 * a compiler warning when we tried to assign it to an object of type
 * UpdateFlags.
 *
 * @ref UpdateFlags
 */
inline
UpdateFlags
operator & (UpdateFlags f1, UpdateFlags f2)
{
  return static_cast<UpdateFlags> (
    static_cast<unsigned int> (f1) &
    static_cast<unsigned int> (f2));
}


/**
 * Global operator which clears all the bits in the first argument if
 * they are not also set in the second argument.
 *
 * @ref UpdateFlags
 */
inline
UpdateFlags &
operator &= (UpdateFlags &f1, UpdateFlags f2)
{
  f1 = f1 & f2;
  return f1;
}








#endif /* UPDATE_FLAGS_HH_ */
