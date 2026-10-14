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
static const std::string TAG_FULLASSEMBLY_CONST = "FullAssembly_const";
static const std::string TAG_FULLASSEMBLY_MODEL = "FullAssembly_model";
static const std::string TAG_FULLASSEMBLY_FE    = "FullAssembly_fe";
static const std::string TAG_COMPUTELOCAL_CONST = "ComputeLocal_const";
static const std::string TAG_COMPUTELOCAL_MODEL = "ComputeLocal_model";
static const std::string TAG_COMPUTELOCAL_FE    = "ComputeLocal_fe";
static const std::string TAG_EVALFIELDS_CONST   = "EvalFields_const";
static const std::string TAG_EVALFIELDS_MODEL   = "EvalFields_model";
static const std::string TAG_EVALFIELDS_FE      = "EvalFields_fe";


#endif /* BENCH_FIELDALGO_INIT_HH_ */
