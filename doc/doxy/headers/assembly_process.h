/**
 * @defgroup assembly_process  Assembly Process 
 *
 * Evaluation of Fields is performed during assemblies at different evaluation points. Assembly process allows summary 
 * evaluation of elements on one patch. The calculation is performed subsequently on all these elements.
 * 
 * 
 * \section problems Solved problems
 *    
 * FieldFE use FEValues in implementation of the `value` method, this involves creation of the FEValues object. This object 
 * should be created only once ideally.
 * 
 * FieldFormula can be evaluated efficiently for the vectors of points. So we need to evaluate lot of points together. That 
 * is possible once the points are from elements of the same region.
 * 
 * We need to achieve effective use of memory:
 *   - efficient ordering of the mesh nodes and elements
 *   - same evaluations of different field types during the assembly process
 *   - use one memory place involved during assembly on the single element
 * Enable advantage of the vector operations.          
 *
 * 
 * \section design Design ideas and constraints
 *
 * The weak formulation is based on the evaluation of the integrals of four types:
 *   - *bulk integral* - decomposed to simplix elements
 *   - *edge integral* - decomposed to edges, edge is set of colocated faces/sides
 *   - *boundary integral* - decomposed to: boundary element of dimension D + faces of bulk elements of dimension D+1 
 *   - *coupling integral* - decomposed to: bulk element of dimension D + faces of elements of dimension D+1 
 * The integration over elementary domain is done in terms of numerical quadrature, i.e. weighted average of the values 
 * in quadrature points. Single logical quadrature point may be distributed to more then one bulk element. 
 *
 * All necessary quadrature points must be distributed to the reference elements during setup of the assmebly. FieldFE is 
 * then able to create and precompute the FEValue object for all these local points. The values on the reference element 
 * are computed once. Mapping (of vectors and tensors) to the actual element can be done only for the active points on current 
 * patch of elements.
 *
 * For the efficient evaluation of the FieldFormula fields we need to evaluate about 128 points at once. Efficient parsing 
 * and formula evaluation library is necessary. We use a specially developed parser for this purpose: *BParser*.
 * 
 * In order to allow FieldFormula to depend on other fields we need that other fields are also evaluated into caches. This is 
 * a bit unfortunate for the FieldConstant and possibly cache inefficient in the case of FieldFE. The FieldFE caching is 
 * inefficient only in the case of using same FE space for some field as well as for the equation itself, that is in the case 
 * of nonlinear equations. Even in such case, we at most compute local base functions twice, once when prefetching the field 
 * value cache and second when assembling on individual elements. If the recalculation takes longer then L3 access we can cache 
 * all FEValue data for all patch quad points.
 * 
 * Example of CPU specification: Intel Core i7: L1d 32kB, L2 256kB, L3 8MB. 
 * [I7 cache specification](https://www.edn.com/design/systems-design/4399725/Memory-Hierarchy-Design---Part-6--The-Intel-Core-i7): 
 * block size 64 bytes, L1d 8way, L2 8way, L3 16way shared for cores  (What if every core works on separate memory chunk,
 * 16 way can be limitting factor?)
 * Considering about 30 scalar fields, 128 points per field and 8 bytes per double the field cache is about 30kB, 
 * that fits (with other data as Dofs, FEValues) well within L2 cache, even if we consider using just half of 
 * it to allow prefetching. Conclusion: size of the field cache that enables amortization of formula evaluation 
 * is about the size suitable for the CPU cache optimization of the assembly process.     
 * 
 * We need to join evaluation of FieldFormula for individual integrals so there will be more complex logic in the selection of field 
 * algorithm since only subset of the fields is involved in every integral.
 * 
 * Gauss quadratures imply no shared quadrature points between bulk and face quadratures. There can be some reuse for boundary and 
 * coupling integrals, but these are relatively rare.
 * 
 * Face integrals should deal with various side permutations.
 * 
 * We should use existing field variables to referencing their cached values.
 * 
 * Although we will have a common set of local evaluation points on the reference element we may allow to use only subset of them on 
 * the actual patch of elements.
 * 
 * The mapping of the pair (element, evaluation point) to the active point on the current patch must be same for all fields invloved 
 * in the assembly. We need a class for this synchronisation.
 * 
 * Ideas to exploit SIMD:
 *   - Once the fields are precomputed the assmbly of local matrix/vector can use SIMD operations as the forms are evaluated uing a same  
 *     formula for every (i-shape, j-shape, q-point). 
 *   - We can possibly even assembly more elements at once. Considering the local matrices/vectors be placed in a continuous block of memory.
 *   - Possible problem is that as we iterate over the shape functions the field values remains constant so we may need to load the field values 
 *     into constant small vectors. 
 *   - FieldFormula - in order to use SIMD we need to organize the patch q-points into the SIMD blocks (4-AVX2, 8-AVX512, seems relevant also 
 *     on GPUs). First atempt will just put elements of the same region into a continuous block, further we may modify patches to optimize number 
 *     of the regions per patch.
 *   - FieldValues - the vectors and tensors has to be transformed be a mapping matrix. One matrix-vector multiplication takes 9 products and 6 
 *     additions these operations should be vectorized over several vectors. Transformation matrices are usually constant on the element so we 
 *     need to repeat them. 
 *   - FieldFE - see FieldValues, further value in a point is linear combination of the same number of the shape functions (for the same dimension 
 *     and order) Vectorisation over points of the elements of the same dimension and order is possible.               
 * 
 * 
 * 
 * @ingroup sssembly
 *
 */
