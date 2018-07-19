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
 * @file    generic_interpolator.hh
 * @brief
 */

#ifndef GENERIC_INTERPOLATOR_HH_
#define GENERIC_INTERPOLATOR_HH_

template <int spacedim, class Value> class FieldFE;


/**
 * Class providing different methods of interpolation between input and output FieldFE.
 *
 */
template <int spacedim, class Value>
class GenericInterpolator
{
public:
	/// Constructor
	GenericInterpolator();

	/// Interpolate data of \p field_in to \p field_out. Use meshes of both fields as source and target mesh.
	void interpolate(FieldFE<spacedim, Value> & field_out, FieldFE<spacedim, Value> & field_in); // TODO maybe add optional parameters
};


#endif /* GENERIC_INTERPOLATOR_HH_ */
