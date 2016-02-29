/*
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 */


/**
 * @file main_doc.hh
 * @mainpage
 *
 *
 * <h1> Flow123d </h1>
 *
 * <h2>Purpose </h2>
 * Flow123d is a simulator of underground water flow and transport processes.
 * It is aimed in particular for large scale simulations and includes models for
 * water flow in fully saturated porous medium, transport of several chemical substances, sorption of
 * chemical substances into the rock matrix, and chemical reactions.
 *
 * <h2>Features </h2>
 * <b> Complex domains </b> The main
 * feature that distinguish it from other similar software is its capability to represent tiny dislocations
 * in large scale domains as 2D and 1D objects and allows interaction of solution on domains of various dimensions.
 *
 * <b> Mixed-hybrid discretization of Darcy flow. </b> The water flow is driven by Darcy flow, this is basically steady state elliptic problem, but
 *  optionally one can consider compressibility of the water and surrounding rock. In the later case we solve a time dependent parabolic problem.
 *  For discretization of the elliptic or parabolic problem we use mixed-hybrid scheme of zero order using lowest order Raviart-Thomas base functions for
 *  discretization of the velocity field and piecewise constant base functions for the pressure and its traces on the interior mesh edges. Main advantage of the
 *  MH scheme is good approximation of the velocity field which is later used in the transport model. For the parabolic case we have implemented a lumping technique
 *  in order to prevent possible oscillations due to violation of the discrete maximum principle.
 *
 * <b> Advection </b> is modeled only by convection (diffusion/dispersion is in developement). We use simple Finite volume scheme with upwind and backward Euler
 *      for time discretization. Unfortunately for highly heterogeneous velocity field which is natural for preferential fracture flow, we have to satisfy CFL condition
 *      which bounds the time step to very small values. Advection module can also compute sorption and dual porosity model (substance exchange between mobile and
 *      immobile pores.
 *
 * <b> Reactions. </b> There are two modules for simulation of chemical reactions. SEMCHEM module can solve nonlinear differential equations rising form
 *      general multicompoent reactiong systems. On the other side this is very costly and is suitable only for small meshes. The second chemical module
 *      is fast but can cope only with linear reactions, i.e. decays.
 *
 * <b> Paralellism. </b> Both the water flow solver and transport solver can run in parallel on distributed memory systems. We use essentially PETSc and MPI libraries.
 *
 * <h2> Main program modules </h2>
 *
 * <b> @ref input_mod "Input" </b> module defines possible formats and structure of input files and their readers.
 *
 * <b> @ref mesh_mod Mesh </b> module contains discretization of a multidimensional computational domain and
 *  geometrical coincidence of mesh elements.
 *
 * <b> @ref flow_mod "Darcy flow" </b> module with mixed-hybrid solver of linear flow equation.
 *
 * <b> @ref transport_mod "Transport" </b> module with model of chemical substances transport.
 *
 * <b> @ref reactions_mod "Reactions" </b> module.
 *
 * @ref Authors "List of Contributors"
 *
 * $LastChangedDate$
 *
 */

/**
 *
 *
 * @page Authors
 *
 * Jan Březina - coordinator, parallelism, schur complements
 *
 * Otto Severýn - original multidimensional flow
 *
 * Milan Hokr - density driven flow
 *
 * Jiří Kopal - transport
 *
 * Jiří Hnídek - GMASH and VTK output classes, infrastructure
 *
 * Jiří Jeníček - JSON reader
 *
 * Lukáš Zedek - SEMCHEM interface, linear reactions
 *
 * Jakub Šístek - two level domain decomposition methods without overlap
 *
 * Dalibor Frydrich
 *
 * Jan Lisal
 *
 * Tomáš Bambuch - profiler class
 *
 * Michal Nekvasil - automatic builds and tests
 */

/**
 * @defgroup system_mod     System module
 * System module contains general support classes for: debugging, error handling, profiling. There are also @ref Vector and @ref VectorId classes with
 * their iterators.
 *
 * @defgroup input_mod      Input
 * This module should contain most of input classes in particular readers for particular file formats. In particular JSON reader, YAML reader, classes
 * allowing definition and access to nodes of input source tree (IST). Module allows definition of input attributes, which are optional for every node
 * of IST. Base attributes are defined in @ref Input::Type::InputAttributes.
 * Defined in Input namespace.
 *
 * @defgroup io_mod         Output
 * This module should contain most of output classes in particular writers for particular file formats. In particular output into GMSH and VTK data formats.
 *
 * @defgroup la_mod         Linear Algebra
 * This module should contain various classes for linear algebra calculations. For small vectors and matrices we would like to use Armadillo library,
 * but meanwhile we have such functionality in @ref math_fce.cc Class @ref LinSys is meant as C++ wrapper for PETSC and possibly for
 * Trilinos. Class @ref SchurComplement provides parallel computation of Schur complements using PETSC library.
 *
 * @defgroup mesh_mod       Mesh
 * This module should contain classes to maintain and access multidimensional mesh with information about coincidence between elements of the meshes.
 * More general we can think about several meshes (even with same dimension) with information about coincidence.
 *
 * @defgroup transport_mod  Advection
 * This module is for advection model. Currently we have only Finite volume implementation without diffusion/dispersion.
 *
 * @defgroup flow_mod       Darcy flow
 * This module contains Mixed-Hybrid and Lumped mixed-hybird discretization of Darcy flow equation for steady and unsteady case. It also contains particular
 * postprocessing functionality as interpolation into continuous finite element space.
 *
 * @defgroup reactions_mod   Chemical Reactions and Decays
 *
 *
 */

#ifndef MAIN_DOC_HH_
#define MAIN_DOC_HH_


#endif /* MAIN_DOC_HH_ */
