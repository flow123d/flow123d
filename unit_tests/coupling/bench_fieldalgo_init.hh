#ifndef BENCH_FIELDALGO_INIT_HH_
#define BENCH_FIELDALGO_INIT_HH_

/**
 * Holds helper classes that allow to dispatch type of FialdAlgo initializing
 * in benchmark tests (FieldConst, FieldModel, FieldFE).
 */


#include <memory>
#include "la/vector_mpi.hh"
#include "fem/dofhandler.hh"


class InitFieldConst {
private:
	/// Forbidden constructor
	InitFieldConst() {}
};

class InitFieldModel {
private:
    /// Forbidden constructor
    InitFieldModel() {}
};

class InitFieldFE {
public:
    static VectorMPI create_data_vec(shared_ptr<DOFHandlerMultiDim> dh, uint seed) {
        VectorMPI data_vec = dh->create_vector();
    	for (uint i=0; i<data_vec.size(); ++i) {
    		data_vec.set( i, 0.05 * ((i+seed)%19 + 1) );
    	}
    	return data_vec;
    }

private:
    /// Forbidden constructor
    InitFieldFE() {}
};

/**
 * Definitions of TAG names
 */
enum AsmTags
{
    fullassembly_const,
    fullassembly_model,
    fullassembly_fe,
    computelocal_const,
    computelocal_model,
    computelocal_fe,
    evalfields_const,
    evalfields_model,
    evalfields_fe
};
static constexpr const char * tag_name(AsmTags tag) {
    constexpr std::array<const char *, 9> tag_names = {
            "FullAssembly_const", "FullAssembly_model", "FullAssembly_fe",
            "ComputeLocal_const", "ComputeLocal_model", "ComputeLocal_fe",
            "EvalFields_const", "EvalFields_model", "EvalFields_fe"
    };
    return tag_names[tag];
}


#endif /* BENCH_FIELDALGO_INIT_HH_ */
