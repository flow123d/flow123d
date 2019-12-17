# Design of field evaluation mechanism

**Problems**
- FieldFE use FEValues in implementation of the `value` method, this involves creation of the FEValues object in every call and mapping points from local to global coordinates before and then back in the method. Moreover the `value` method is called several times per single element.
- FieldFormula can be evaluated efficiently for the vectors of points (see [muparserx](https://beltoforion.de/article.php?a=muparserx&s=idFeatures#idFeatures) and [exprtk](https://github.com/ArashPartow/exprtk) and the parser [comparison](https://github.com/ArashPartow/math-parser-benchmark-project). So we need to evaluate lot of points together. That is possible once the points are from elements of the same region.
- Current assembly process have probably poor memory locality due to:
  1. inefficient ordering of the mesh nodes and elements
  2. various field evaluation during the assembly process
  3. to many memory places involved during assembly on the single element
  
  

  

**Design ideas and constraints**
- The weak formulation is based on the evaluation of the integrals of four types: 
    - *bulk integral* - decomposed to simplix elements
    - *face integral* - decomposed to faces/sides of elements
    - *edge integral* - decomposed to edges, edge is set of colocated faces/sides
    - *coupling integral* - decomposed to: bulk element of dimension D + faces of elements of dimension D+1
  The integration over elementary domain is done in terms of numerical quadrature, i.e. weighted average of the values in quadrature points. Single logical quadrature point may be distributed to more then one bulk element.
  
- All necessary quadrature points must be distributed to the reference elements during setup of the assmebly.
  FieldFE is then able to create and precompute the FEValue object for all these local points. 
  The values on the reference element are computed once. Mapping (of vectors and tensors) to the actual element
  can be done only for the active points on current patch of elements.
- For the efficient evaluation of the FieldFormula fields we need to evaluate about 128 points at once. 
  Efficient parsing and formula evaluation library is necessary. ExprTK is fast but lacks several important features:
  - SIMD evaluation
  - support of sparse evaluation on the vector data (includes size change)
  - matrix multiplication (matrix-matrix, matrix-vector, vector-vector = dot product)
  Seems that no parser have this kind of functionality. Ultimately, we have to modify ExprTK.
- In order to allow FieldFormula to depend on other fields we need that other fields are also evaluated into caches. 
  This is a bit unfortunate for the FieldConstant and possibly cache inefficient in the case of FieldFE. 
  The FieldFE caching is inefficient only in the case of using same FE space for some field as well as 
  for the equation itself, that is in the case of nonlinear equations. Even in such case, we at most 
  compute local base functions twice, once when prefetching the field value cache and second when 
  assembling on individual elements. If the recalculation takes longer then L3 access we can cache 
  all FEValue data for all patch quad points.

  Example of CPU specification: Intel Core i7: L1d 32kB, L2 256kB, L3 8MB. 
  [I7 cache specification](https://www.edn.com/design/systems-design/4399725/Memory-Hierarchy-Design---Part-6--The-Intel-Core-i7): 
  block size 64 bytes, L1d 8way, L2 8way, L3 16way shared for cores  (What if every core works on separate memory chunk,
  16 way can be limitting factor?)
  Considering about 30 scalar fields, 128 points per field and 8 bytes per double the field cache is about 30kB, 
  that fits (with other data as Dofs, FEValues) well within L2 cache, even if we consider using just half of 
  it to allow prefetching. Conclusion: size of the field cache that enables amortization of formula evaluation 
  is about the size suitable for the CPU cache optimization of the assembly process. 
- We need to join evaluation 
  of FieldFormula for individual integrals so there will be more complex logic in the selection of field algorithm
  since only subset of the fields is involved in every integral. 
-  ~If these quadrature schemes share the local points the values may be evaluated just once.~
  Gauss quadratures imply no shared quadrature points between bulk and face quadratures. There can be 
  some reuse for boundary and coupling integrals, but these are relatively rare. 
- Face integrals should deal with various side permutations.
- We should use existing field variables to referencing their cached values.
- Although we will have a common set of local evaluation points on the reference element we may allow 
  to use only subset of them on the actual patch of elements.   
- The mapping of the pair (element, evaluation point) to the active point on the current patch must 
  be same for all fields invloved in the assembly. We need a class for this synchronisation.


 
# Design (TODO: Update rest according to the Cache structure)

## Caching structures
**EvalPoints** 
Summary all quadrature points on the reference elements. Necessary for the FieldFE precalculation.
- Internaly single table per dimension, but we need to distribute quadrature points accross dimension so there has 
  to be one shared object. However quadratures are added per dimension in the Assembly<dim> classes.
- Every reference element (dimension) have own list of local points. Only resulting 'Integral' objects know ranges of 
  their points.

**Integral**
Every kind of integral is associated with an "assembly object"
- BulkIntegral - on Cell
- EdgeIntegral - on edge.side_range
- CouplingIntegral - on bulk.ngh_range
- BoundaryIntegral - on boundary_side_range (TBD)

- In theory we can combine integrals to have at most single instance of each type per single assembly operation. 
- For every dimension we create them by a method of the EvalPoints, passing suitable set of quadratures.
- Functions provided by Integrals:
  - Mark active eval points on the actual patch. Assembly algorithm iterates over "elementary integration domains" 
    in the order of their Hilbert index. The integral has to iterate over its quadrature points. The point  and distributing the points on them until we reach the number of active points per patch.
- Same method Integral::points(assembly_object) is used later on to iterate over quadrature points.
- In principle there should be one kind of the quadrature point accessor for every kind of the integral.

Operations provided by the quadrature points:
	- all points - retrieve value of a field in the point, e.g. `pressure(bulk_point)`		
	- SidePoint - permute(CellSide) -> SidePoint; allows iteration over all pairs of sides of an edge, in fact allows
	  to form integrals mixing points of any two faces of same dimension, even non-colocated.
	- CouplingSidePoint - bulk() -> BulkPoint
Possible interface to shape functions:
	`FieldFE pressure_fe ...`
	`element_fields.normal(side_point)` - outer normal at a side point
	`element_fields.coord(point)` - absolute coordinates of the point
	`element_fields.JxW(point)` - Jacobina times quadrature weight of the point
    `pressure_fe.base(i, point)` - value of the i-th shape function in the point
	`pressure_fe.grad(i, point)` - value of the i-th base function in the point

**ElementCacheMap**
This holds mesh elements associated with the actual patch and for every pair (mesh element, eval_point) gives the active point index, these have to form continuous sequence. The points of the "assemble objects" forming the patch are added to the map. There is independent map for every integral but these can be in a single class. 
There should be no active points reuse ase the active points of the two "assemble objects" of the same kind are disjoint (can change if we decide to change elemetaty "assembly objects", e.g. singel NGHSide instead of the range).
ElementCacheMap is shared by dimensions !! It should be part of the EqData.

**FieldValueCache**
- One per field, holds only values for selected integrals.
- One value is a scalar, vector, or tensor. But stored through Array in singel chunk on memory.

1. Every Field<spacedim, Value> has its cache for its values. 
    The logical cache elements have type given by the Value template parameter however as this is inefficient for the Armadillo objects we allocate plain memory and construct the value objects referencing to the allocated memory.
2. Every Field<> has an array of three instances of the FieldValueCache one for each dimension. We should try to make FieldValueCache not templated by the dimension to avoid virtual calls during the read access to the cache.
3. FieldValueCache is logically a table of N element slots composed of M values for the M local points on the reference element. The implementation uses a plain memory or a fixed size array.

## Cache operations
**Overview**
Usage of the field caches consists of:
1. Merging more quadratures into a single set of the local evaluation points (class `EvalPoints`).
2. Create a `FieldSet` one for every quadrature of the fields involved in that integral. 
3. Initialize the fields in the integral's field set: Allocate the cache space in fields and mark which quadrature the field use. This is done through the call of `FieldSet::cache_allocate<dim>(EvalSubset, Mapping<dim>)`
4. Composing the assambly patch, element cache prefetching. This is organized by 
`ElementCacheMap`.
5. Assembly: Iterate over Integrals and get cached values from the fields.


### 1. Initialization - evaluation points

**EvalPoints**
The class to store set of common local points. 
Two operations:
1. Add a point set, return their indices in the table. ~Eliminate duplicate points 
    with prescribed tolerance.~ Contrary, we want to keep the points from the single quadrature in a single continuous block.
    In fact, we probably need to use two distinguised methods:
	`BulkIntegral &&EvalPoints::add_bulk(Quadrature bulk_quadrature)`
	`BoundaryIntegral &&add_boundary(Quadrature side_quadrature)`
	`EdgeIntegral &&add_edge(Quadrature side_quadrature)`
	`CouplingIntegral &&add_coupling(Quadrature bulk_quadrature)`
	
2. Get local point coordinates for given index.

Common storage of quadrature points is not necessary anymore, but we keep it in order to keep indexing of the points the same accross the fields (which can use ony certain subblocks of the points). We can decrease the memory footprint (and thus the cache footprint) at the expense of slower access to the cached values (one more indirect memory access).

**EvalSubset**
The object containing array of indices into local point set, there is only single array for the 
bulk quadrature, but (n_sides x n_permutations) for a side quadrature. We shall try to keep this 
common for the bulk and side quadratures and non-templated. It keeps pointer to the EvalPoints, in particular to know total number of local points.


**Example**
```
class Assembly<dim> {
	EvalPoints ep;
	EvalSubset face_eval;
	
	Assembly(EqData &eq) {
	    // Distribute quadrature points to the reference elements through EvalPoints
	    EvalPoints &ep = eq.ep;
	    ep.add_bulk(Gauss(dim, order));
		ep.add_edge(Gauss(dim-1, order));
		ep.add_boundary(Gauss(dim-1, order));
		ep.add_ngh(Gauss(dim, order));
			
			Gauss(dim, order)
		eEvalSubset stiffness_eval;
	...
	
	
	Mapping map;
	map.setup_mapping_fields(this->all_fields);
	ep.set_mapping(map);
	....
	this->mass_eval = this->ep.add_bulk();
	this->face_eval = this->ep.add_side(Gauss(dim-1, order));
        ...
```

### 2. Integration field sets
This use existing field sets to simplify group operations on those. 
```
    // still in class Assembly<dim> definition
    bulk_fields = eqdata.subset({'cross_section', 'porosity'})
    edge_fields = eqdata.subset({'cross_section', 'flux'})
	...
	
	// cache allocation
	bulk_fields.cache_allocate(eval_points.bulk_points<dim>().)
```

### 3. Initialization - cache allocation
`FieldSet::cache_allocate` just calls the same method for every `FieldCommon` in the set.

```
FieldCommon::cache_allocate(EvalSubset sub_quad)
```
This is a virtual method with implementation in Field<...> which allocates `FieldValueCache<Value>(sub_quad.eval_points())` during the first call and mark the used local points `FieldValueCache<Value>::mark_used(sub_quad)`. 

The Field<> have array of these classes, one instance for every dimension so the `FieldValueCache` 
must not be templated (by dimension). The only field algorithm that needs absolute coordinates is FieldFormula. We will pass the whole EqData fieldset to it in order to allow more complex dependencies. One of the fields will be 'Coordinate' or 'Position' field which just compute absolute coordinates of the local points (using the mapping). 
`FieldValueCache` has a vector with start and size of active blocks of the quadrature points.

TODO: Different assembly loops may need different evaluation subsets, so either we need mean to reallocate caches before every assembly loop or just update only selected part of the field value cache. Current structure of FieldValueCache allows to have allocated individual EvalSubset tables separately and update only selected subsets. Moreover the other subsets do not block the CPU cache (example of associativity, seems that AMD processors have 4way L1, 16way L2 and 64way L3 cache; try to get similar info for Intel) anyway we shuld mak an analysis how many memory blocks we need in any computation loop. Innermost loops should idealy access at most 4 memory blocks (4 pointers). 

### FieldValueCache
In principle this is just a table of items of type Value with dimensions: n_cached_elements x n_evaluation_points. Other properties to keep:
- reference to the EvalPoints instance (therefore it should not be dim templated either)
- dimension (probably just for checks)
- mask which points (of the EvalPoints) are used in passed `EvalSubset` objects, only values in these points has to be evaluated.

**Cache structure**
- FieldFormula can update in single step all values for common subset and common region
- so we want to have such values in continuous blocks in the cache
- Cache is organized into blocks of the size (n_points x n_max_elements)  for every subset, where n_points is number of points in the subset and n_max_elements is max number of elements in the single patch.
- Every cache block is indexed as `[i_element][i_point]`, since within the block we always update all points, but only elements on the same region, and moreover we reuse already cached elements.


**Example**
```
    // still in class Assembly<dim> definition
   this->mass_fields.cache_allocate(this->mass_eval);
   this->side_fields.cache_allocate(this->side_eval);
} // end of Assembly<dim>
```
TODO:
- FieldValueCache allocates its table at the first call, but the mask for active local points is added from more calls


### 4. ElementCacheMap
New algo:
1. Set whole table [i_cache_element][i_eval_point] to -1.
2. Iterate through the elementary integrals, add them to the actual patch until we reach given number of points.
   The elements to which we distribute the points are added to the element_map, mapping global element ID to its cache index.
   Add active points of single elementary integral, i.e. set [i_cache_element][i_eval_point] to the n_active_points. 
3. Pass through the every field. For single field:
   1. Pass through the patch elements:
      ```
	  for (el : el_cache_map.patch_elements())
	  		// Some algorithms may update the field value cache right here
			if (field.region_fields[el.region_id].mark_cache_points(el_cache_map.points(el))) {
				postponed_fields_set.add(field.region_fields[el.region_id]);
	  ```
   2. Update marked points, release marks.
   ```
   for (f in postponed_fields_set)
   		f.fill_cache(el_cache);
   ```		   	
      		
			

This class synchronize the cached elements between (all) fields of single equation. It provides mapping from elements to the cache index and list of elements to cache. This have overloaded evaluation operator, which returns the index in the cache for the given element (or element index). The last cache line is overwritten if the index is not in the cache.
The implementation use: table cache_idx -> el_idx, hash mapping el_idx -> cache_idx, list of cache lines that schould be updated.
```
void Assembly::mass_assembly(DHCellAccessor cell) {
    // fill element cache index in the DHCellAccessor
    DHCellAccessor el_cache = this->element_cache_map(cell);    

    // More elements can be cached in the generic assembly loop as 
    // the cache prefetching.
    
```
Proposed interface:

`int ElementCacheMap::operator ()(cell.element())`
Maps element to its cache line, throuws if it is not in the cache.

`int ElementCacheMap::add(cell.element())`
Add an element to the cache (if not presented), replace the oldes cache line.

`vector<int> added_elements;`
Public vector of cache lines to update. Indices appended by the `add` method. FieldSet::cache_update clean the list
after all fields are updated.


#### Cache Update
```
    // Call cache_update for the fields in the field sets
    // This can also be done in the generic loop.
    this->eq_data.cache_update(this->element_cache_map);    
    
    // Possible update of the base function values.
    presssure_field_fe.fe_values.update(cell);
```
Two major algorithms are in use:
- FieldFE - evaluates base func values in all quadrature points (done once per assembly),  dot product with DOFs, optionaly multiplied by the Mapping matrix (important optimization for vector fields and derivatives, must have support in FEValues)
- FieldFormula - evaluates all elements in the patch (same region), in all point from single continuous block od quad points


### 5 Assembly, cache read
```
    /*
    // Asumme following types:
    EvalSubset this->mass_eval;
    EvalSubset this->side_eval;
    EvalSubset this->ngh_side_eval;
    */
    
    ...
    DHCellAccessor cache_cell = this->element_cache_map(cell);    
    // Bulk integral, no sides, no permutations.
    for(BulkPoint q_point: this->mass_eval.points(cache_cell)) {
        // Extracting the cached values.
        double cs = cross_section(q_point);        
        
	// Following would be nice to have. Not clear how to 
        // deal with more then single element as fe_values have its own cache that has to be updated.
        auto base_fn_grad = presssure_field_fe.base_value(q_point);
	loc_matrix += outer_product((cs * base_fn_grad),  base_fn_grad)
    }   

    // Side integrals.
    // FieldFE<..> conc;
    for (DHCellSide side : cache_cell.side_range()) {        
	for(DHCellSide el_ngh_side : side.edge_sides()) {            
   	    // vector of local side quadrature points in the correct side permutation
	    Range<SidePoint> side_points = this->side_eval.points(side)
	    for (SidePoint p : side_points) {
	    	ngh_p = p.permute(el_ngh_side);
	        loc_mat += cross_section(p) * sigma(p) * 
		    (conc.base_value(p) * velocity(p) 
		    + conc.base_value(ngh_p) * velocity(ngh_p)) * p.normal() / 2;
            }
       }
    }
    
    // Dimension coupling
    // TODO: update
    for (DHNeighbSide ngh_side cache_cell.neighb_sides()) {        
        // vector of local side quadrature points in the correct side permutation
        Range<BulkPoint> side_points = this->ngh_side_eval.points(ngh_side);
        for (auto p : side_points) {
            side_avg += cross_section(el_ngh_p) * sigma(el_ngh_p) * 
                ( velocity(side_p) + velocity(el_ngh_p)) * side_p.normal();
        }
    }
}
```

In order to use FieldValueCache in consistent way, we can precompute: quadrature points coordinates, normals and JxW into suitable FieldValueCache and possibly wrap it into Fields in order to use them in FieldFormula. However direct access should be provided through BulkPoint and SidePoint. **precise**
We want cache_allocate and cache_update to update the mapping caches however we have to pass the mapping somehow into XYPoint and to do this on single place we want to pass it to EvalPoints.

**Interface for the bulk integrals**

Bulk integrals are evaluated only on the single element, so we can assume, that all values are cached. No need for 
evaluationg just some points etc.
`Field<..>::operator() (BulkPoint)`
use dimension and cache element index that are part of the BulkPoint structure to access the cached value.
Similarly we can exted caching in FeValues to more elements and provide acces to them:
`FieldFE<..>::base_value(BulkPoint, component)`

`FieldFE<..>::base_grad(BulkPoint, component)`

**Interface for the face integrals*

We need to access values on two elements of two matching sides. All sides of all elements are necessary
so we want to cache all side points in the same way as the values in bulk points. We can access matching side
but we also need its proper permutation, to this end we provide SidePoint and its mapping to the connected side of the other element.

`Field<..>::operator() (SidePoint)`
`FieldFE<..>::base_value(SidePoint, component)`
`FieldFE<..>::base_grad(SidePoint, component)`

**Interface for dimension coupling integrals**
In general the quadrature can be different then the quadrature used on faces, or the could be no integration over the faces as in the case of P1 method. No problem for the lower dim element as we have to evaluate fields at all bulk points on al these elements. Only matching evaluation points on the connected side of a bulk element are necessary. Moreover we are not able to change mask of evaluation point accoring to the elements. 

Two possible solutions:
- local point sets varies with elements
- evaluation of noncached values, but still want to have the FieldValues in FieldFE for their points



TODO: 
- No sense to have more then single bulk element in bulk integrals.
- In Side integrals need: normals, more then single element (elements of same or different dimension).
- How to match q poiunts on side of higher dim element and bulk of lower dim elemnet. Not clear how it works in current implementation. Is the side permutation compatible with connected fracture elements?
- How to cache values in quadrature points on faces of the higher dim elements efficiently? We want to precompute only values on a single side. It could be fine in case of DG if we use same quadratures for side terms as well as for the coupling temrs. But for general quadratures there could be problem. It seems the either we do not cache in this case or need more general cache.
- How to apply only FeValues to FE living on edges without FeSideValues.




**Further thoughts**
- Extension to interdependent fields: If a field depends on the other field, it recursively informs the other field about quadrature.

- We introduce fields for absolute cooridinates X, Y, Z  as well as for the depth, this is related to the generalization of the FieldFormula, that can use also other fields in the formulas.



