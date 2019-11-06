#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest_mpi.hh>

#include "system/global_defs.h"



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
            y_v[i] = (double)i/100;
            z_v[i] = (10 - (double)(i%10)) * 1.1;
        }
    }

	void profiler_output() {
		static ofstream os( FilePath("parsers_test.log", FilePath::output_file) );
		Profiler::instance()->output(MPI_COMM_WORLD, os);
		os << endl;
	}

    double fast_compute_constant() {
    	double sum = 0.0;
        START_TIMER("1_fast_compute");
        for (int j=0; j<nLoops; ++j) {
            for (int i=0; i<nBulkSize; ++i) {
                result_v[i] = 0.5;
            }
        }
        END_TIMER("1_fast_compute");
        for (int i=0; i<nBulkSize; ++i) {
        	sum += result_v[i];
        }

    	return sum;
    }

    double fast_compute_simple() {
    	double sum = 0.0;
        START_TIMER("1_fast_compute");
        for (int j=0; j<nLoops; ++j) {
            for (int i=0; i<nBulkSize; ++i) {
                result_v[i] = x_v[i] + y_v[i] + z_v[i];
            }
        }
        END_TIMER("1_fast_compute");
        for (int i=0; i<nBulkSize; ++i) {
            sum += result_v[i];
        }

    	return sum;
    }

    double fast_compute_complex() {
    	double sum = 0.0;
    	double pi = 3.141592653589793238462643;
        START_TIMER("1_fast_compute");
        for (int j=0; j<nLoops; ++j) {
            for (int i=0; i<nBulkSize; ++i) {
                result_v[i] = 2*x_v[i] + y_v[i]*3 + x_v[i]*(z_v[i]-y_v[i]) + 2*pi*z_v[i];
            }
        }
        END_TIMER("1_fast_compute");
        for (int i=0; i<nBulkSize; ++i) {
            sum += result_v[i];
        }

    	return sum;
    }

    double fast_compute_power() {
    	double sum = 0.0;
    	double pi = 3.141592653589793238462643;
        START_TIMER("1_fast_compute");
        for (int j=0; j<nLoops; ++j) {
            for (int i=0; i<nBulkSize; ++i) {
                result_v[i] = 2*x_v[i] + pow( y_v[i], 3.0 ) + x_v[i]*(z_v[i]-y_v[i]) + 2*pi*z_v[i];
            }
        }
        END_TIMER("1_fast_compute");
        for (int i=0; i<nBulkSize; ++i) {
            sum += result_v[i];
        }

    	return sum;
    }

    double exprtk_parse_substitution(std::string expression_string) {
        typedef exprtk::symbol_table<double> symbol_table_t;
        typedef exprtk::expression<double>   expression_t;
        typedef exprtk::parser<double>       parser_t;

        double res = 0.0;
        double x, y, z;

        symbol_table_t symbol_table;
        symbol_table.add_variable("x",x);
        symbol_table.add_variable("y",y);
        symbol_table.add_variable("z",z);
        symbol_table.add_constants();

        expression_t expression;
        expression.register_symbol_table(symbol_table);

        parser_t parser;
        parser.compile(expression_string,expression);

        START_TIMER("2_exprtk_parse_substitution");
        for (int j=0; j<nLoops; ++j)
        	for (int i=0; i<nBulkSize; ++i)
            {
        		x = x_v[i]; y = y_v[i]; z = z_v[i];
        		result_v[i] = expression.value();
            }
        END_TIMER("2_exprtk_parse_substitution");
    	for (int i=0; i<nBulkSize; ++i)
        {
            res += result_v[i];
        }

    	return res;
    }

    /*double exprtk_parse_vector(std::string expr) {
        typedef exprtk::symbol_table<double> symbol_table_t;
        typedef exprtk::expression<double>     expression_t;
        typedef exprtk::parser<double>             parser_t;

        std::string expression_string =
                      " for (var i := 0; i < min(x[],y[],z[],r[]); i += 1) { r[i] := " + expr + "  }";

        symbol_table_t symbol_table;
        symbol_table.add_vector("x",x_v);
        symbol_table.add_vector("y",y_v);
        symbol_table.add_vector("z",z_v);
        symbol_table.add_vector("r",result_v);
        symbol_table.add_constants();

        expression_t expression;
        expression.register_symbol_table(symbol_table);

        parser_t parser;
        parser.compile(expression_string,expression);

        double res = 0.0;
        START_TIMER("exprtk_parse_vector");
        for (int j=0; j<nLoops; ++j) {
            expression.value();
        }
        END_TIMER("exprtk_parse_vector");
        for (int i=0; i<nBulkSize; ++i) res += result_v[i];

    	return res;
    }*/

    /*double exprtk_parse_vector_with_view(std::string expr) {
        typedef exprtk::symbol_table<double> symbol_table_t;
        typedef exprtk::expression<double>     expression_t;
        typedef exprtk::parser<double>             parser_t;

        std::string expression_string =
                      " for (var i := 0; i < min(x[],y[],z[],r[]); i += 1) { r[i] := " + expr + "  }";

        exprtk::vector_view<double> x_view = exprtk::make_vector_view(x_v,x_v.size());
        exprtk::vector_view<double> y_view = exprtk::make_vector_view(y_v,y_v.size());
        exprtk::vector_view<double> z_view = exprtk::make_vector_view(z_v,z_v.size());
        exprtk::vector_view<double> r_view = exprtk::make_vector_view(result_v,result_v.size());

        symbol_table_t symbol_table;
        symbol_table.add_vector("x",x_view);
        symbol_table.add_vector("y",y_view);
        symbol_table.add_vector("z",z_view);
        symbol_table.add_vector("r",r_view);
        symbol_table.add_constants();

        expression_t expression;
        expression.register_symbol_table(symbol_table);

        parser_t parser;
        parser.compile(expression_string,expression);

        double res = 0.0;
        START_TIMER("5_exprtk_parse_vector_with_view");
        for (int j=0; j<nLoops; ++j) {
            expression.value();
        }
        END_TIMER("5_exprtk_parse_vector_with_view");
        for (int i=0; i<nBulkSize; ++i) res += result_v[i];

    	return res;
    }*/

    /*double exprtk_parse_vector_view_rebase(std::string expression_string) {
        typedef exprtk::symbol_table<double> symbol_table_t;
        typedef exprtk::expression<double>     expression_t;
        typedef exprtk::parser<double>             parser_t;

        double sum = 0.0;
        std::vector<double> v(3, 0.0);
        exprtk::vector_view<double> view = exprtk::make_vector_view(v,v.size());

        symbol_table_t symbol_table;
        symbol_table.add_vector("v",view);
        symbol_table.add_constants();

        expression_t expression;
        expression.register_symbol_table(symbol_table);

        parser_t parser;
        parser.compile(expression_string,expression);

        int i;
        START_TIMER("6_exprtk_parse_vector_view_rebase");
        for (int j=0; j<nLoops; ++j) {
            i=0;
        	for (auto& new_vec : vv)
            {
                view.rebase(new_vec.data()); // update vector
                result_v[i] += expression.value();
            }
        	i++;
        }
        END_TIMER("6_exprtk_parse_vector_view_rebase");
        for (int i=0; i<nBulkSize; ++i) sum += result_v[i];

        return sum;
    }*/

    double exprtk_parse_vector_fast(std::string expr) {
        typedef exprtk::symbol_table<double> symbol_table_t;
        typedef exprtk::expression<double>     expression_t;
        typedef exprtk::parser<double>             parser_t;

        std::string expression_string = " r := " + expr + " ";

        exprtk::vector_view<double> x_view = exprtk::make_vector_view(x_v,x_v.size());
        exprtk::vector_view<double> y_view = exprtk::make_vector_view(y_v,y_v.size());
        exprtk::vector_view<double> z_view = exprtk::make_vector_view(z_v,z_v.size());
        exprtk::vector_view<double> r_view = exprtk::make_vector_view(result_v,result_v.size());

        double sum = 0.0;

        symbol_table_t symbol_table;
        symbol_table.add_vector("x",x_view);
        symbol_table.add_vector("y",y_view);
        symbol_table.add_vector("z",z_view);
        symbol_table.add_vector("r",r_view);
        symbol_table.add_constants();

        expression_t expression;
        expression.register_symbol_table(symbol_table);

        parser_t parser;
        parser.compile(expression_string,expression);

        expression.value();

        START_TIMER("3_exprtk_parse_vector_fast");
        for (int j=0; j<nLoops; ++j) {
            expression.value();
        }
        END_TIMER("3_exprtk_parse_vector_fast");
        for (int i=0; i<nBulkSize; ++i) sum += result_v[i];

        return sum;
    }

    void run_tests(int vec_size) {
        create_data_vectors(vec_size);
        std::cout << "Tests is running with " << vec_size << " points ..." << std::endl;

        START_TIMER("A_constant_expresions");
        this->fast_compute_constant();
        this->exprtk_parse_substitution(constantLine);
        this->exprtk_parse_vector_fast(constantLine);
        END_TIMER("A_constant_expresions");

        START_TIMER("B_simple_expresions");
        this->fast_compute_simple();
        this->exprtk_parse_substitution(simpleLine);
        this->exprtk_parse_vector_fast(simpleLine);
        END_TIMER("B_simple_expresions");

        START_TIMER("C_complex_expresions");
        this->fast_compute_complex();
        this->exprtk_parse_substitution(complexLine);
        this->exprtk_parse_vector_fast(complexLine);
        END_TIMER("C_complex_expresions");

        START_TIMER("D_power_expresions");
        this->fast_compute_power();
        this->exprtk_parse_substitution(powerLine);
        this->exprtk_parse_vector_fast(powerLine);
        END_TIMER("D_power_expresions");
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
};


TEST(Parser, all) {
    ParserHandler pHandler;

    START_TIMER("test_16_points");
    pHandler.run_tests(16);
    END_TIMER("test_16_points");

    START_TIMER("test_32_points");
    pHandler.run_tests(32);
    END_TIMER("test_32_points");

    START_TIMER("test_64_points");
    pHandler.run_tests(64);
    END_TIMER("test_64_points");

    START_TIMER("test_128_points");
    pHandler.run_tests(128);
    END_TIMER("test_128_points");

    START_TIMER("test_256_points");
    pHandler.run_tests(256);
    END_TIMER("test_256_points");

    START_TIMER("test_512_points");
    pHandler.run_tests(512);
    END_TIMER("test_512_points");

    START_TIMER("test_1024_points");
    pHandler.run_tests(1024);
    END_TIMER("test_1024_points");

    START_TIMER("test_2048_points");
    pHandler.run_tests(2048);
    END_TIMER("test_2048_points");

    pHandler.profiler_output();
}

/* ***************************************************************************************************************************** */
/*TEST(Parser, exprtk_vector) {
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double>     expression_t;
    typedef exprtk::parser<double>             parser_t;

    std::string expression_string = "2*v[0] + v[1]^3 + v[0]*(v[2]-v[1]) + 2*pi*v[2]";

    double sum = 0.0;
    std::vector<double> v(3, 0.0);
    std::vector<std::vector<double>> vv = { {1.1, 1.0, 2.2}, {2.2, 2.0, 3.3}, {3.3, 3.0, 4.4}, {4.4, 4.0, 5.5}, {5.5, 5.0, 6.6} };
    exprtk::vector_view<double> view = exprtk::make_vector_view(v,v.size());

    symbol_table_t symbol_table;
    symbol_table.add_vector("v",view);
    symbol_table.add_constants();

    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;
    parser.compile(expression_string,expression);

    for (auto& new_vec : vv)
    {
       view.rebase(new_vec.data()); // update vector
       sum = expression.value();
       printf("result: %19.15f\n", sum);
    }

    //for (unsigned int i=0; i<5; ++i) printf("result: %19.15f\n", res[i]);
}*/
