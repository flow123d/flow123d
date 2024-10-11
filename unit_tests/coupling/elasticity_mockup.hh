#ifndef ELASTICITY_MOCKUP_HH_
#define ELASTICITY_MOCKUP_HH_


#include <mesh_constructor.hh>
#include "arma_expect.hh"
#include <rev_num.h>

#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "fields/field_value_cache.hh"
#include "fields/field_values.hh"
#include "fields/field_set.hh"
#include "fields/field_fe.hh"
#include "fields/generic_field.hh"
#include "fields/multi_field.hh"
#include "fields/bc_multi_field.hh"
#include "fields/equation_output.hh"
#include "fields/field_model.hh"
#include "fields/field_constant.hh"
#include "coupling/equation.hh"
#include "tools/unit_si.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"
#include "fem/fe_p.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "la/linsys.hh"
#include "la/linsys_PETSC.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"
#include "system/sys_profiler.hh"


class ElasticityMockupTest : public testing::Test {
public:
	ElasticityMockupTest()
    {
		string root_dir=string(UNIT_TESTS_BIN_DIR) + "/coupling";
		string build = string(__DATE__) + ", " + string(__TIME__)
	            + " flags: (unknown compiler flags)";

        FilePath::set_io_dirs(".",root_dir,"",".");
        Profiler::instance();
        Profiler::instance()->set_program_info("Flow123d",
                string(FLOW123D_VERSION_NAME_), string(FLOW123D_GIT_BRANCH_), string(FLOW123D_GIT_REVISION_), build);
        Profiler::set_memory_monitoring(false, false);
    }

    ~ElasticityMockupTest() {}

    /// Perform profiler output.
    void profiler_output(std::string file_name) {
		FilePath fp(file_name + "_profiler.json", FilePath::output_file);
		Profiler::instance()->output(MPI_COMM_WORLD, fp.filename());
	}
};

#endif /* ELASTICITY_MOCKUP_HH_ */
