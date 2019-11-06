#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest_mpi.hh>

#include "system/global_defs.h"

#ifdef FLOW123D_RUN_UNIT_BENCHMARKS


#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <iostream>
#include <limits>
#include <ios>
#include <iomanip>
#include <numeric>

#include "system/sys_profiler.hh"
#include "system/file_path.hh"

#include "fields/exprtk.hpp"

using namespace std;


class ParserHandler {
public:
    ParserHandler() {
	    Profiler::initialize();
	    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

	    constantLine = "0.5";
        simpleLine = "x + y + z";
        complexLine = "2*x + y*3 + x*(z-y) + 2*pi*z";
        powerLine = "2*x + y^3 + x*(z-y) + 2*pi*z";

        funcPlus = "y+0.5";
        funcPower = "y^3";
        funcAbs = "abs(z)";
        funcExp = "exp(y)";
        funcLog = "log(y)";
        funcSin = "sin(pi*y)";
        funcAsin = "asin(z)";
        funcTernary = "z>0 ? x : y";
    	funcMax = "max(x,y,z)";
    }

    ~ParserHandler() {
	    Profiler::uninitialize();
    }

    void create_data_vectors(int vec_size) {
        nBulkSize = vec_size;
        nLoops = 5000 * (2048 / vec_size);

        x_v.resize(nBulkSize);
        y_v.resize(nBulkSize);
        z_v.resize(nBulkSize);
        result_v.resize(nBulkSize);

        for (int i=0; i<nBulkSize; ++i) {
            x_v[i] = i;
            y_v[i] = (double)i/100 + 0.01;
            z_v[i] = (4.95 - (double)(i%10) * 1.1) * 0.2; // values oscillate in <-1; 1>
        }
    }

	void profiler_output(std::string file_name) {
		static ofstream os( FilePath(file_name, FilePath::output_file) );
		Profiler::instance()->output(MPI_COMM_WORLD, os);
		os << endl;
	}

    double cpp_compute_constant() {
    	double sum = 0.0;
        START_TIMER("cpp_compute");
        for (int j=0; j<nLoops; ++j) {
            for (int i=0; i<nBulkSize; ++i) {
                result_v[i] = 0.5;
            }
        }
        END_TIMER("cpp_compute");
        for (int i=0; i<nBulkSize; ++i) {
        	sum += result_v[i];
        }

    	return sum;
    }

    double cpp_compute_simple() {
    	double sum = 0.0;
        START_TIMER("cpp_compute");
        for (int j=0; j<nLoops; ++j) {
            for (int i=0; i<nBulkSize; ++i) {
                result_v[i] = x_v[i] + y_v[i] + z_v[i];
            }
        }
        END_TIMER("cpp_compute");
        for (int i=0; i<nBulkSize; ++i) {
            sum += result_v[i];
        }

    	return sum;
    }

    double cpp_compute_complex() {
    	double sum = 0.0;
    	double pi = 3.141592653589793238462643;
        START_TIMER("cpp_compute");
        for (int j=0; j<nLoops; ++j) {
            for (int i=0; i<nBulkSize; ++i) {
                result_v[i] = 2*x_v[i] + y_v[i]*3 + x_v[i]*(z_v[i]-y_v[i]) + 2*pi*z_v[i];
            }
        }
        END_TIMER("cpp_compute");
        for (int i=0; i<nBulkSize; ++i) {
            sum += result_v[i];
        }

    	return sum;
    }

    double cpp_compute_complex_power() {
    	double sum = 0.0;
    	double pi = 3.141592653589793238462643;
        START_TIMER("cpp_compute");
        for (int j=0; j<nLoops; ++j) {
            for (int i=0; i<nBulkSize; ++i) {
                result_v[i] = 2*x_v[i] + pow( y_v[i], 3.0 ) + x_v[i]*(z_v[i]-y_v[i]) + 2*pi*z_v[i];
            }
        }
        END_TIMER("cpp_compute");
        for (int i=0; i<nBulkSize; ++i) {
            sum += result_v[i];
        }

    	return sum;
    }

    double exprtk_parse_vector_fast(std::string expr) {
        typedef exprtk::symbol_table<double> symbol_table_t;
        typedef exprtk::expression<double>     expression_t;
        typedef exprtk::parser<double>             parser_t;

        std::string expression_string = " result_vec := " + expr + " ";

        exprtk::vector_view<double> x_view = exprtk::make_vector_view(x_v,x_v.size());
        exprtk::vector_view<double> y_view = exprtk::make_vector_view(y_v,y_v.size());
        exprtk::vector_view<double> z_view = exprtk::make_vector_view(z_v,z_v.size());
        exprtk::vector_view<double> r_view = exprtk::make_vector_view(result_v,result_v.size());

        double sum = 0.0;

        symbol_table_t symbol_table;
        symbol_table.add_vector("x",x_view);
        symbol_table.add_vector("y",y_view);
        symbol_table.add_vector("z",z_view);
        symbol_table.add_vector("result_vec",r_view);
        symbol_table.add_constants();

        expression_t expression;
        expression.register_symbol_table(symbol_table);

        parser_t parser;
        parser.compile(expression_string,expression);

        expression.value();

        START_TIMER("exprtk_parse_vector_fast");
        for (int j=0; j<nLoops; ++j) {
            expression.value();
        }
        END_TIMER("exprtk_parse_vector_fast");
        for (int i=0; i<nBulkSize; ++i) sum += result_v[i];

        return sum;
    }

    void run_expression_tests(int vec_size) {
        create_data_vectors(vec_size);
        std::cout << "Tests is running with " << vec_size << " points ..." << std::endl;

        START_TIMER("A_constant_expresions");
        this->cpp_compute_constant();
        this->exprtk_parse_vector_fast(constantLine);
        END_TIMER("A_constant_expresions");

        START_TIMER("B_simple_expresions");
        this->cpp_compute_simple();
        this->exprtk_parse_vector_fast(simpleLine);
        END_TIMER("B_simple_expresions");

        START_TIMER("C_complex_expresions");
        this->cpp_compute_complex();
        this->exprtk_parse_vector_fast(complexLine);
        END_TIMER("C_complex_expresions");

        START_TIMER("D_power_expresions");
        this->cpp_compute_complex_power();
        this->exprtk_parse_vector_fast(powerLine);
        END_TIMER("D_power_expresions");
        std::cout << " ... OK" << std::endl;
    }

    void run_function_tests(int vec_size) {
        create_data_vectors(vec_size);
        std::cout << "Tests is running with " << vec_size << " points ..." << std::endl;

        START_TIMER("A_plus_function");
        this->exprtk_parse_vector_fast(funcPlus);
        END_TIMER("A_plus_function");

        START_TIMER("B_power_function");
        this->exprtk_parse_vector_fast(funcPower);
        END_TIMER("B_power_function");

        START_TIMER("C_abs_function");
        this->exprtk_parse_vector_fast(funcAbs);
        END_TIMER("C_abs_function");

        START_TIMER("D_exp_function");
        this->exprtk_parse_vector_fast(funcExp);
        END_TIMER("D_exp_function");

        START_TIMER("E_log_function");
        this->exprtk_parse_vector_fast(funcLog);
        END_TIMER("E_log_function");

        START_TIMER("F_sin_function");
        this->exprtk_parse_vector_fast(funcSin);
        END_TIMER("F_sin_function");

        START_TIMER("G_asin_function");
        this->exprtk_parse_vector_fast(funcAsin);
        END_TIMER("G_asin_function");

        START_TIMER("H_ternary_function");
        this->exprtk_parse_vector_fast(funcTernary);
        END_TIMER("H_ternary_function");

        START_TIMER("I_max_function");
        this->exprtk_parse_vector_fast(funcMax);
        END_TIMER("I_max_function");
        std::cout << " ... OK" << std::endl;
    }

    // data members
    int nLoops;
    int nBulkSize;

    std::vector<double> x_v;
    std::vector<double> y_v;
    std::vector<double> z_v;
    std::vector<double> result_v;

    std::string constantLine;
    std::string simpleLine;
    std::string complexLine;
    std::string powerLine;

    std::string funcPlus;
    std::string funcPower;
    std::string funcAbs;
    std::string funcExp;
    std::string funcLog;
    std::string funcSin;
    std::string funcAsin;
    std::string funcTernary;
	std::string funcMax;
};


/**
 * Speed test of base expressions parsing in exprtk:
 *  - constant
 *  - simple
 *  - complex
 *  - complex with power function
 *
 * All tests are compared with evaluation of same expression specified by raw C++ code.
 */
TEST(Parser, expressions) {
	// test of base expressions: constant, simple, complex and complex with power function
    ParserHandler pHandler;

    START_TIMER("test_16_points");
    pHandler.run_expression_tests(16);
    END_TIMER("test_16_points");

    START_TIMER("test_32_points");
    pHandler.run_expression_tests(32);
    END_TIMER("test_32_points");

    START_TIMER("test_64_points");
    pHandler.run_expression_tests(64);
    END_TIMER("test_64_points");

    START_TIMER("test_128_points");
    pHandler.run_expression_tests(128);
    END_TIMER("test_128_points");

    START_TIMER("test_256_points");
    pHandler.run_expression_tests(256);
    END_TIMER("test_256_points");

    START_TIMER("test_512_points");
    pHandler.run_expression_tests(512);
    END_TIMER("test_512_points");

    START_TIMER("test_1024_points");
    pHandler.run_expression_tests(1024);
    END_TIMER("test_1024_points");

    START_TIMER("test_2048_points");
    pHandler.run_expression_tests(2048);
    END_TIMER("test_2048_points");

    pHandler.profiler_output("exprtk_expressions.yaml");
}


/**
 * Speed test of selected functions parsing in exprtk:
 *  - plus (comparative operation for all other)
 *  - power
 *  - abs
 *  - exp
 *  - log
 *  - sin
 *  - asin
 *  - ternary operator
 *  - max
 */
TEST(Parser, functions) {
    ParserHandler pHandler;

    START_TIMER("test_128_points");
    pHandler.run_function_tests(128);
    END_TIMER("test_128_points");

    START_TIMER("test_256_points");
    pHandler.run_function_tests(256);
    END_TIMER("test_256_points");

    START_TIMER("test_512_points");
    pHandler.run_function_tests(512);
    END_TIMER("test_512_points");

    START_TIMER("test_1024_points");
    pHandler.run_function_tests(1024);
    END_TIMER("test_1024_points");

    pHandler.profiler_output("exprtk_functions.yaml");
}

#endif // FLOW123D_RUN_UNIT_BENCHMARKS
