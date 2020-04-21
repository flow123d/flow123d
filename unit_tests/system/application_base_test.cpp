
#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>

#include "system/application_base.hh"
#include "system/system.hh"


class ApplicationBaseTest : public testing::Test, public ApplicationBase {
public:
	ApplicationBaseTest() : testing::Test(), ApplicationBase() {}
protected:
	void run() {
		xprintf(Err, "testing error...\n");
	}

    void parse_cmd_line(const int, char **) override 
    {}
	
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
	EXPECT_THROW_WHAT( {run();}, ExcXprintfMsg, "testing error...");
    EXPECT_THROW_WHAT( {seg_fault();}, ExcSignal, "Signal 11" );
}
