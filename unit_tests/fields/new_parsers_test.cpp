#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest_mpi.hh>

#include "system/global_defs.h"


#include "_parsers/muParserTest.h"

/** \brief This macro will enable mathematical constants like M_PI. */
#define _USE_MATH_DEFINES

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

#include "_parsers/muParser.h"
#include "_parsers/exprtk.hpp"

using namespace std;
using namespace mu;


class ParserHandler {
public:
    ParserHandler() {
	    Profiler::initialize();
	    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

	    constantLine = "0.5";
        simpleLine = "x + y + z";
        complexLine = "2*x + y^3 + x*(z-y) + 2*_pi*z";
        complexLine_pi = "2*x + y^3 + x*(z-y) + 2*pi*z";
        constantLine_i = "0.5";
        simpleLine_i = "x[i] + y[i] + z[i]";
        complexLine_i = "2*x[i] + y[i]^3 + x[i]*(z[i]-y[i]) + 2*pi*z[i]";
        constantLine_v = "0.5";
        simpleLine_v = "v[0] + v[1] + v[2]";
        complexLine_v = "2*v[0] + v[1]^3 + v[0]*(v[2]-v[1]) + 2*pi*v[2]";
    }

    ~ParserHandler() {
	    Profiler::uninitialize();

	    dealloc_data_vectors();
    }

    void dealloc_data_vectors() {
        delete [] x;
        delete [] y;
        delete [] z;
        delete [] result;
    }

    void create_data_vectors(int vec_size) {
	    nBulkSize = vec_size;
	    x = new value_type[nBulkSize];
    	y = new value_type[nBulkSize];
    	z = new value_type[nBulkSize];
    	result = new value_type[nBulkSize];

    	x_v.resize(nBulkSize);
    	y_v.resize(nBulkSize);
    	z_v.resize(nBulkSize);
    	vv.resize(nBulkSize);
    	result_v.resize(nBulkSize);

        for (int i=0; i<nBulkSize; ++i) {
            x[i] = i;
            y[i] = (value_type)i/100;
            z[i] = (10 - (value_type)(i%10)) * 1.1;
            x_v[i] = x[i];
            y_v[i] = y[i];
            z_v[i] = z[i];
            vv[i].resize(3);
            vv[i][0] = x[i]; vv[i][1] = y[i]; vv[i][2] = z[i];
        }
    }

	void profiler_output() {
		static ofstream os( FilePath("parsers_test.log", FilePath::output_file) );
		Profiler::instance()->output(MPI_COMM_WORLD, os);
		os << endl;
	}

    value_type fast_compute_constant() {
        value_type sum = 0.0;
        START_TIMER("1_fast_compute");
        for (int j=0; j<nLoops; ++j) {
            for (int i=0; i<nBulkSize; ++i) {
                result[i] = 0.5;
            }
        }
        END_TIMER("1_fast_compute");
        for (int i=0; i<nBulkSize; ++i) {
        	sum += result[i];
        }

    	return sum;
    }

    value_type fast_compute_simple() {
        value_type sum = 0.0;
        START_TIMER("1_fast_compute");
        for (int j=0; j<nLoops; ++j) {
            for (int i=0; i<nBulkSize; ++i) {
                result[i] = x[i] + y[i] + z[i];
            }
        }
        END_TIMER("1_fast_compute");
        for (int i=0; i<nBulkSize; ++i) {
            sum += result[i];
        }

    	return sum;
    }

    value_type fast_compute_complex() {
        value_type sum = 0.0;
        value_type pi = 3.141592653589793238462643;
        START_TIMER("1_fast_compute");
        for (int j=0; j<nLoops; ++j) {
            for (int i=0; i<nBulkSize; ++i) {
                result[i] = 2*x[i] + y[i]*y[i]*y[i] + x[i]*(z[i]-y[i]) + 2*pi*z[i];
            }
        }
        END_TIMER("1_fast_compute");
        for (int i=0; i<nBulkSize; ++i) {
            sum += result[i];
        }

    	return sum;
    }

    value_type mu_parse_expresion(string_type sLine) {
        value_type sum = 0.0;
        mu::Parser  parser;
   	    parser.DefineVar(_T("x"), x);
        parser.DefineVar(_T("y"), y);
        parser.DefineVar(_T("z"), z);
        parser.SetExpr(sLine);
        START_TIMER("2_mu_parser");
        for (int j=0; j<nLoops; ++j) {
            parser.Eval(result, nBulkSize);
        }
        END_TIMER("2_mu_parser");
        for (int i=0; i<nBulkSize; ++i) sum += result[i];

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

        START_TIMER("3_exprtk_parse_substitution");
        for (int j=0; j<nLoops; ++j)
        	for (int i=0; i<nBulkSize; ++i)
            {
        		x = x_v[i]; y = y_v[i]; z = z_v[i];
        		result[i] = expression.value();
            }
        END_TIMER("3_exprtk_parse_substitution");
    	for (int i=0; i<nBulkSize; ++i)
        {
            res += result[i];
        }

    	return res;
    }

    double exprtk_parse_vector(std::string expr) {
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
        START_TIMER("4_exprtk_parse_vector");
        for (int j=0; j<nLoops; ++j) {
            expression.value();
        }
        END_TIMER("4_exprtk_parse_vector");
        for (int i=0; i<nBulkSize; ++i) res += result_v[i];

    	return res;
    }


    double exprtk_parse_vector_with_view(std::string expr) {
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
    }

    double exprtk_parse_vector_view_rebase(std::string expression_string) {
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
    }

    double exprtk_parse_vector_fast(std::string expr) {
        typedef exprtk::symbol_table<double> symbol_table_t;
        typedef exprtk::expression<double>     expression_t;
        typedef exprtk::parser<double>             parser_t;

        std::string expression_string = " r := " + expr + " ";

        double sum = 0.0;

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

        expression.value();

        START_TIMER("7_exprtk_parse_vector_fast");
        for (int j=0; j<nLoops; ++j) {
            expression.value();
        }
        END_TIMER("7_exprtk_parse_vector_fast");
        for (int i=0; i<nBulkSize; ++i) sum += result_v[i];

        return sum;
    }

    void run_tests(int vec_size) {
        create_data_vectors(vec_size);
        std::cout << "Tests is running with " << vec_size << " points ..." << std::endl;

        START_TIMER("A_constant_expresions");
        this->fast_compute_constant();
        this->mu_parse_expresion(constantLine);
        this->exprtk_parse_substitution(constantLine_i);
        //this->exprtk_parse_vector(constantLine_i);
        //this->exprtk_parse_vector_with_view(constantLine_i);
        //this->exprtk_parse_vector_view_rebase(constantLine_v);
        this->exprtk_parse_vector_fast(constantLine);
        END_TIMER("A_constant_expresions");

        START_TIMER("B_simple_expresions");
        this->fast_compute_simple();
        this->mu_parse_expresion(simpleLine);
        this->exprtk_parse_substitution( (std::string)simpleLine );
        //this->exprtk_parse_vector( simpleLine_i );
        //this->exprtk_parse_vector_with_view( simpleLine_i );
        //this->exprtk_parse_vector_view_rebase(simpleLine_v);
        this->exprtk_parse_vector_fast(simpleLine);
        END_TIMER("B_simple_expresions");

        START_TIMER("C_complex_expresions");
        this->fast_compute_complex();
        this->mu_parse_expresion(complexLine);
        this->exprtk_parse_substitution(complexLine_pi);
        //this->exprtk_parse_vector(complexLine_i);
        //this->exprtk_parse_vector_with_view(complexLine_i);
        //this->exprtk_parse_vector_view_rebase(complexLine_v);
        this->exprtk_parse_vector_fast(complexLine_pi);
        END_TIMER("C_complex_expresions");
        std::cout << " ... OK" << std::endl;
    }

    // data members
    static const int nLoops = 10000;

    int nBulkSize;

    value_type *x;
    value_type *y;
    value_type *z;
    value_type *result;

    std::vector<double> x_v;
    std::vector<double> y_v;
    std::vector<double> z_v;
    std::vector<std::vector<double>> vv;
    std::vector<double> result_v;

    string_type constantLine;
    string_type simpleLine;
    string_type complexLine;
    std::string complexLine_pi;
    std::string constantLine_i;
    std::string simpleLine_i;
    std::string complexLine_i;
    std::string constantLine_v;
    std::string simpleLine_v;
    std::string complexLine_v;
};


TEST(Parser, all) {
    ParserHandler pHandler;

    START_TIMER("test_0250_points");
    pHandler.run_tests(250);
    END_TIMER("test_0250_points");
    pHandler.dealloc_data_vectors();

    START_TIMER("test_0500_points");
    pHandler.run_tests(500);
    END_TIMER("test_0500_points");
    pHandler.dealloc_data_vectors();

    START_TIMER("test_1000_points");
    pHandler.run_tests(1000);
    END_TIMER("test_1000_points");
    pHandler.dealloc_data_vectors();

    START_TIMER("test_2000_points");
    pHandler.run_tests(2000);
    END_TIMER("test_2000_points");

    pHandler.profiler_output();
} // */

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

/*TEST(Parser, exprtk_vector) {
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double>     expression_t;
    typedef exprtk::parser<double>             parser_t;

    std::string expression_string =
                  " var pi := 3.141592;                                       "
                  " for (var i := 0; i < min(x[],y[],z[],r[]); i += 1)        "
                  " {                                                         "
                  "   r[i] := 2*x[i] + y[i]^3 + x[i]*(z[i]-y[i]) + 2*pi*z[i]; "  //"2*x[i] + y[i]^3 + x[i]*(z[i]-y[i]) + 2*pi*z[i]";
                  " }                                                         ";

    std::vector<double> x = { 1.1, 2.2, 3.3, 4.4, 5.5 };
    std::vector<double> y = { 1.0, 2.0, 3.0, 4.0, 5.0 };
    std::vector<double> z = { 2.2, 3.3, 4.4, 5.5, 6.6 };
    std::vector<double> res = { 0.0, 0.0, 0.0, 0.0, 0.0 };

    //exprtk::vector_view<double> x_view = exprtk::make_vector_view(x,x.size());
    //exprtk::vector_view<double> y_view = exprtk::make_vector_view(y,y.size());
    //exprtk::vector_view<double> z_view = exprtk::make_vector_view(z,z.size());
    //exprtk::vector_view<double> r_view = exprtk::make_vector_view(res,res.size());

    symbol_table_t symbol_table;
    symbol_table.add_vector("x",x);
    symbol_table.add_vector("y",y);
    symbol_table.add_vector("z",z);
    symbol_table.add_vector("r",res);
    symbol_table.add_constants();

    expression_t expression;
    expression.register_symbol_table(symbol_table);

    parser_t parser;
    parser.compile(expression_string,expression);

    expression.value();
    //printf("result: %19.15f\n", ret);

    for (unsigned int i=0; i<5; ++i) printf("result: %19.15f\n", res[i]);
} // */
