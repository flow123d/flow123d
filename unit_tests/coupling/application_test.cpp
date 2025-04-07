
#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>

#include "system/system.hh"
#include "coupling/application.hh"
#include "moc_mpi_bench_test_error.hh"

class ApplicationBaseTest : public testing::Test, public Application {
public:
	ApplicationBaseTest() : testing::Test(), Application() {}
protected:
	void run() {
		ASSERT_PERMANENT(false).error("testing error...\n");
	}

	void seg_fault() {
      // Attempt to read from unallocated memory.
      petsc_initialize(0, nullptr);
      int *i = nullptr;
      printf("%d", i[0]);
    }

    void SetUp() override
    {}
    void TearDown() override
    {}
};


TEST_F(ApplicationBaseTest, Exceptions) {
	EXPECT_THROW_WHAT( {run();}, feal::Exc_assert, "testing error...");
    EXPECT_THROW_WHAT( {seg_fault();}, ExcSignal, "Signal 11" );
}
