#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <vector>

#include "system/global_defs.h"


#ifdef FLOW123D_RUN_UNIT_BENCHMARKS

#include "system/armor.hh"
#include "system/sys_profiler.hh"

static const uint N = 100000;
static const uint REPEAT = 1000;

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


TEST(Armor_speed, compare_speed) {
    
    Profiler::instance();
    
    std::vector<arma::vec3> armaArray1(N);
    for (auto & item : armaArray1) {
        item[0] = fRand(-100,100);
        item[1] = fRand(-100,100);
        item[2] = fRand(-100,100);
    }
    
    std::vector<arma::vec3> armaArray2(N);
    for (auto & item : armaArray2) {
        item[0] = fRand(-100,100);
        item[1] = fRand(-100,100);
        item[2] = fRand(-100,100);
    }
    
    double armaDotResults[N];
    
    arma::vec3 armaSum{0, 0, 0};
    
    arma::Mat<double>::fixed<3, 3> armaMat{1, 2, 3, 4, 5, 6, 7, 8, 9};
    
    std::vector<Armor::Mat<double, 3, 1>> armorArray1(N);
    for (auto & item : armorArray1) {
        item[0] = fRand(-100,100);
        item[1] = fRand(-100,100);
        item[2] = fRand(-100,100);
    }
    
    std::vector<Armor::Mat<double, 3, 1>> armorArray2(N);
    for (auto & item : armorArray2) {
        item[0] = fRand(-100,100);
        item[1] = fRand(-100,100);
        item[2] = fRand(-100,100);
    }
    
    double armorDotResults[N];
    
    Armor::Mat<double, 3, 1> armorSum{0, 0, 0};
    
    Armor::Mat<double, 3, 3> armorMat{1, 2, 3, 4, 5, 6, 7, 8, 9};
    
    {
        START_TIMER("armadillo_sum");
        for (uint j{0}; j < REPEAT; ++j) {
            for (uint i{0}; i < N; ++i) {
                armaSum = armaSum + armaArray1[i];
            }
        }
        END_TIMER("armadillo_sum");
    }
    {
        START_TIMER("armor_sum");
        for (uint j{0}; j < REPEAT; ++j) {
            for (uint i{0}; i < N; ++i) {
                armorSum = armorSum + armorArray1[i];
            }
        }
        END_TIMER("armor_sum");
    }
    {
        START_TIMER("armadillo_dot");
        for (uint j{0}; j < REPEAT; ++j) {
            for (uint i{0}; i < N; ++i) {
                armaDotResults[i] = dot(armaArray1[i], armaArray2[i]);
            }
        }
        END_TIMER("armadillo_dot");
    }
    {
        START_TIMER("armor_dot");
        for (uint j{0}; j < REPEAT; ++j) {
            for (uint i{0}; i < N; ++i) {
                armorDotResults[i] = dot(armorArray1[i], armorArray2[i]);
            }
        }
        END_TIMER("armor_dot");
    }
    {
        START_TIMER("arma_multiplication");
        for (uint j{0}; j < REPEAT; ++j) {
            for (uint i{0}; i < N; ++i) {
                armaArray2[i] = armaMat * armaArray1[i];
            }
        }
        END_TIMER("arma_multiplication");
    }
    {
        START_TIMER("armor_multiplication");
        for (uint j{0}; j < REPEAT; ++j) {
            for (uint i{0}; i < N; ++i) {
                armorArray2[i] = armorMat * armorArray1[i];
            }
        }
        END_TIMER("armor_multiplication");
    }
    /*
     * Results need to be used somehow so the compiler doesn't skip unnecessary calculations
    */
    auto armaForPrinting{armaArray2[0]};
    auto armorForPrinting{armorArray2[0]};
    EXPECT_EQ(0, ((armaForPrinting[0] + armorForPrinting[0] + armaSum[0] + armorSum[0] + armaDotResults[0] + armorDotResults[0]) * 0));
    Profiler::instance()->output(cout);
    Profiler::uninitialize();
}

#endif // FLOW123D_RUN_UNIT_BENCHMARKS
