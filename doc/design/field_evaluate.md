# Design of field evaluation mechanism

**Problems**
- FieldFE use FEValues in implementation of the `value` method, this involves creation of the FEValues object in every call and mapping points from local to global coordinates before and then back in the method. Moreover the `value` method is called several times per single element.
- FieldFormula can be evaluated efficiently for the vectors of points (see [muparserx](https://beltoforion.de/article.php?a=muparserx&s=idFeatures#idFeatures) and [exprtk](https://github.com/ArashPartow/exprtk) and the parser [comparison](https://github.com/ArashPartow/math-parser-benchmark-project). So we need to evaluate lot of points together. That is possible once the points are from elements of the same region.

**Design ideas and constraints**
- Have a sort of field value cache for the local points on the reference element.
- The assembly process may use several different integrals using several quadrature schemes. One bulk and one boundary integral is the minimal choice. ~If these schemes share the same local points the values may be evaluated just once.~
Gauss quadratures imply no shared quadrature points, moreover continuous blocks of points are necessary for efficient FieldFormula implementation.
- Allocation of the cache space and setup of quadratures should be done once per assembly.
- The cache should cover several elements in order to exploit shared evaluation points on the element sides. The cache should be updated on several elements at once to amortize overhead in FieldFormula parser. 
- Boundary quadratures should deal with various side permutations.
- We should use existing field variables to referencing their cached values.
- Although we will have a common set of local evaluation points on the reference element
  different fields may update only a subset of these values according to the integrals in which they appear.
- The mapping of cache slots to actual mesh elements must be same for all fields invloved in the assembly. We need a mean for this synchronisation.


 
# Design 

## Cache structure
1. Every Field<spacedim, Value> has its cache for its values. 
    The logical cache elements have type given by the Value template parameter however as this is inefficient for the Armadillo objects we allocate plain memory and construct the value objects referencing to the allocated memory.
2. Every Field<> has an array of three instances of the FieldValueCache one for each dimension. We should try to make FieldValueCache not templated by the dimension to avoid virtual calls during the read access to the cache.
3. FieldValueCache is logically a table of N element slots composed of M values for the M local points on the reference element. The implementation uses a plain memory or a fixed size array.

## Cache operations
**Overview**
Usage of the field caches consists of:
1. Merging more quadratures into a single set of the local evaluation points (class `EvalPoints`).
    For every merged quadrature we obtain the EvalSubset. 
2. Create a `FieldSet` one for every quadrature of the fields involved in that integral. 
3. Initialize the fields in the integral's field set: Allocate the cache space in fields and mark which quadrature the field use. This is done through the call of `FieldSet::cache_allocate<dim>(EvalSubset, Mapping<dim>)`
4. In the main assembly loop the element cache prefetching can be done. This is organized by 
`ElementCacheMap` which knows which elements are cached.
5. Assembly on a single element:
    1. Update caches of the used fields (can be possibly moved into generic assembly loop).
    2. Map the element (elements) to their cache indices.
    3. Iterate over EvalSubsets and get cached values from the fields.


### 1. Initialization - evaluation points

**EvalPoints**
The class to store set of common local points. 
Two operations:
1. Add a point set, return their indices in the table. ~Eliminate duplicate points 
    with prescribed tolerance.~ Contrary, we want to keep the points from the single quadrature in a single continuous block.
    In fact, we probably need to use two distinguised methods:
	`EvalSubset add_bulk(Quadrature bulk_quadrature)`
	`EvalSubset add_side(Quadrature side_quadrature)`
2. Get point for given index.

Common storage of quadrature points is not necessary anymore, but we keep it in order to keep indexing of the points the same accross the fields (which can use ony certain subblocks of the points). We can decrease the memory footprint (and thus the cache footprint) at the expense of slower access to the cached values (one more indirect memory access).

**EvalSubset**
The object containing array of indices into local point set, there is only single array for the 
bulk quadrature, but (n_sides x n_permutations) for a side quadrature. We shall try to keep this 
common for the bulk and side quadratures and non-templated. It keeps pointer to the EvalPoints, in particular to know total number of local points.


**Example**
```
class Assembly<dim> {
	EvalPoints ep;
	EvalSubset mass_eval;
	EvalSubset face_eval;
	EvalSubset stiffness_eval;
	....
	this->mass_eval = this->ep.add_bulk(Gauss(dim, order));
	this->face_eval = this->ep.add_side(Gauss(dim-1, order));
        ...
```

### 2. Integration field sets
This use existing field sets to simplify group operations on those. 
```
    // still in class Assembly<dim> definition
    this->mass_fields = eqdata.subset({'cross_section', 'porosity'})
    this->face_fields = eqdata.subset({'cross_section', 'flux'})
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

### FieldValueCache
In principle this is just a table of items of type Value with dimensions: n_cached_elements x n_evaluation_points. Other properties to keep:
- reference to the EvalPoints instance (therefore it should not be dim templated either)
- dimension (probably just for checks)
- mask which points (of the EvalPoints) are used in passed `EvalSubset` objects, only values in these points has to be evaluated.


**Example**
```
    // still in class Assembly<dim> definition
   this->mass_fields.cache_allocate(this->mass_eval);
   this->side_fields.cache_allocate(this->side_eval);
} // end of Assembly<dim>
```
TODO:
- FieldValueCache allocates its table at the first call, but the mask for active local points is added from more calls


### 5.2 ElementCacheMap
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


### 5.1 (and 4.) Cache Update
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

### 5.3 Cache read
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



