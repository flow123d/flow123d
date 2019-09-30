# Design of field evaluation mechanism

**Problems**
- FieldFE use FEValues in implementation of the `value` method, this involves creation of the FEValues object in every call and mapping points from local to global coordinates before and then back in the method. Moreover the `value` method is called several times per single element.
- FieldFormula can be evaluated efficiently for the vectors of points (see [muparserx](https://beltoforion.de/article.php?a=muparserx&s=idFeatures#idFeatures) and [exprtk](https://github.com/ArashPartow/exprtk) and the parser [comparison](https://github.com/ArashPartow/math-parser-benchmark-project). So we need to evaluate lot of points together. That is possible once the points are from elements of the same region.

**Design ideas and constraints**
- Have a sort of field value cache for the local points on the reference element.
- The assembly process may use several different integrals using several quadrature schemes. One bulk and one boundary integral is the minimal choice. If these schemes share the same local points the values may be evaluated just once.
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
    The logical cache elements have type given by the Value template parameter however as this is inefficient for the Armadillo objects we allocate plain memory and construct the value object over its patches.
2. Every Field<> has an array of three instances of the FieldValueCache one for each dimension. We should try to make FieldValueCache not templated by the dimension to avoid virtual calls during the read access to the cache.
3. FieldValueCache is logically a table of N element slots composed of M values for the M local points on the reference element. The implementation uses a plain memory or a fixed size array.

## Cache operations
**Overview**
Usage of the field caches consists of:
1. Merging more quadratures into a single set of the local evaluation points (class `EvalPoints`).
    For every merged quadrature we obtain the EvalSubset. 
2. Create a `FieldSet` one for every qudrature of the fields involved in that integral. 
3. Initialize the fields in the integral's field set: Allocate the cache space in fields and mark 
the local points used in integral's quadrature. This is done through the call of `FieldSet::cache_allocate<dim>(EvalSubset, Mapping<dim>)`
4. In the main assembly loop the element cache prefetching can be done. This is organized by 
`ElementCacheMap` this knows which elements are cached.
5. Assembly on a single element:
    1. Update caches of the used fields (can be possibly moved into generic assembly loop).
    2. Map the element (elements) to their cache indices.
    3. Iterate over EvalSubsets and get cached values from the fields.


### 1. Initialization - evaluation points

**EvalPoints**
The class to store set of common local points. 
Two operations:
1. Add a point set, return their indices in the table. Eliminate duplicate points 
    with prescribed tolerance.
    In fact, we probably need to use two distinguised methods:
	`EvalSubset add_subset<dim>(Quadrature<dim> bulk_quadrature)`
	`EvalSubset add_subset<dim>(SideQuadrature<dim> side_quadrature)`
2. Get point for given index.

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
	this->mass_eval = this->ep.add_subset(Gauss<dim>(order));
	this->face_eval = this->ep.add_subset(SideQuadrature(Gauss<dim-1>(order)));
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
FieldCommon::cache_allocate(vector<EvalSubset> sub_quads, Mixed<Mapping> mapping)
```
This is a virtual method with implementation in Field<...> which stores the mapping in the field and calls  `FieldValueCache<Value>(sub_quads)`. The Field<> have array of these classes, one instance for every dimension so the class should not be templated (by dimension)

**FieldValueCache**
In principle this is just a table of items of type Value with dimensions: n_cached_elements x n_evaluation_points. Other properties to keep:
- reference to the EvalPoints instance (therefore it should not be dim templated either)
- dimension (probably just for checks)
- mask which points (of the EvalPoints) are used in passed `EvalSubset` objects, only values in these points has to be evaluated.
- allocation informations are just collected, actual allocation is performed during first update call
  alternatively we can introduce explicit function call to do that.

**Example**
```
    // still in class Assembly<dim> definition
   this->mass_fields.cache_allocate(this->mass_eval, this->mapping);
   this->side_fields.cache_allocate(this->side_eval, this->mapping);
} // end of Assembly<dim>
```
TODO:
- FieldValueCache allocates its table at the first call, but the mask for active local points is added from more calls

### 5.2 Map element to the cache element index
```
void Assembly::mass_assembly(DHCellAccessor cell) {
    // this can possibly be called also in the generic assembly loop as 
    // the cache prefetching
    // assume that ngh_el, is a neigbour element
    
    // map used element into element cache
    el_cache = this->element_cache_map(cell.element());
    ngh_el_cache = this->element_cache_map(ngh_el.element());

```
### 5.1 (and 4.) Cache Update
```
    // just call cache_update for the fields in the field sets
    this->mass_fields.cache_update(this->element_cache_map);
    this->face_fields.cache_update(this->element_cache_map);
```
### 5.3 Cache read
```
    // Bulk integral, no sides, no permutations.
    for(auto q_point: this->mass_eval.points(el_cache)) {
        // Extracting the cached values.
        double cs = cross_section(q_point);
        double ngh_cs = cross_section.value(ngh_el_cache, q_point);
        // This would be nice to have. Not clear how to 
        // deal with more then single element as fe_values have its own cache that has to be updated.
        double base_fn = presssure_field_fe.fe_values.value(q_point);
    }   

    // Side integrals.
    for (auto side : cell.side_range()) {
        double side_avg = 0;
		for(auto el_ngh_side : side.edge_sides()) {            
   	    	// vector of local side quadrature points in the correct side permutation
	        auto side_points = this->side_eval.points(side->side())
		    for (auto p : side_points) {
		        side_avg += cross_section(side, p) * sigma(side, p) * 
				    ( velocity(side, p) + velocity(el_ngh_side, p)) * p.normal();
            }
       }
    }
    
    // Dimension coupling
    for (auto ngh_side cell.neighb_sides()) {        
        // vector of local side quadrature points in the correct side permutation
        auto ngh_side_points = this->ngh_side_eval.points(ngh_side)
        for (auto p : side_points) {
            side_avg += cross_section(side, p) * sigma(side, p) * 
                ( velocity(side, p) + velocity(el_ngh_side, p)) * p.normal();
        }
    }
}
```
- `Field<spacedim, Value>::value(EvalElement eval_el, EvalPoint q_point)`
  returns Value, i.e. some FieldValue<..>. The implementation should do something like:
  ```return this->value_cache[eval_el.dim()].value(eval_el.cache_idx(), q_point.cache_idx())```
  Where `value_cache` is the instance of `FieldValueCache`.
- `this->mass_eval` and `this->side_eval` are instances of `EvalSubset`.


TODO: 
- No sense to have more then single bulk element in bulk integrals.
- In Side integrals need: normals, more then single element (elements of same or different dimension).
- How to match q poiunts on side of higher dim element and bulk of lower dim elemnet. Not clear how it works in current implementation. Is the side permutation compatible with connected fracture elements?
- How to cache values in quadrature points on faces of the higher dim elements efficiently? We want to precompute only values on a single side. It could be fine in case of DG if we use same quadratures for side terms as well as for the coupling temrs. But for general quadratures there could be problem. It seems the either we do not cache in this case or need more general cache.
- How to apply only FeValues to FE living on edges without FeSideValues.




**Further thoughts**
- Extension to interdependent fields: If a field depends on the other field, it recursively informs the other field about quadrature.

- We introduce fields for absolute cooridinates X, Y, Z  as well as for the depth, this is related to the generalization of the FieldFormula, that can use also other fields in the formulas.



