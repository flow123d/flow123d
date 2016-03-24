
#include <flow_gtest.hh>

#include "system/application_base.hh"
#include "system/system.hh"


class ApplicationBaseTest : public testing::Test, public ApplicationBase {
public:
	ApplicationBaseTest() : ApplicationBase(0, NULL), testing::Test() {}
protected:
	void run() {
		xprintf(Err, "testing error...\n");
	}
	
	void seg_fault() {
      // Attempt to read from unallocated memory.
      petsc_initialize(0, nullptr);
      int *i = nullptr;
      printf("%d", i[0]);
    }

    virtual void SetUp() {
    }
    virtual void TearDown() {
    };
};


TEST_F(ApplicationBaseTest, Exceptions) {
	EXPECT_THROW_WHAT( {run();}, ExcXprintfMsg, "testing error...");
    EXPECT_THROW_WHAT( {seg_fault();}, ExcSignal, "Signal 11" );
}