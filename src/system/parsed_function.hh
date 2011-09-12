/*
 * FunctionParser.hh
 *
 *  Created on: Sep 8, 2011
 *      Author: jb
 */

#ifndef FUNCTIONPARSER_HH_
#define FUNCTIONPARSER_HH_

#include "fparser.hh"
#include <string>


class ParsedFunction {
public:
    ParsedFunction() : changed(true) {}

    void set_expression(const std::string &expr) {
        expression=expr;
        changed=true;
    }

    inline void set_constant(std::string& name, double value) {
        ASSERT(parser.AddConstant(name, value), "Illegal constant name in parsed function!\n");
        changed=true;
    }

    inline void set_time(double value) {
        parser.AddConstant("t", value);
        changed=true;
    }

    inline double value(arma::vec3& point)
    {
        if (changed) parse();
        double res=parser.Eval(point.memptr());
        DBGMSG("point: %f %f %f val: %f\n",point.memptr()[0],point.memptr()[1],point.memptr()[2],res);
        return res;
    }
private:
    void parse()
    {
        std::string vars = "x,y,z";
        DBGMSG("%s::%s\n",expression.c_str(),vars.c_str());
        parser.Parse(expression, vars);
        INPUT_CHECK(parser.GetParseErrorType() == FunctionParser::FP_NO_ERROR,"Parsed Function: %s\n",parser.ErrorMsg());
        parser.Optimize();
        changed=false;
    }

    bool changed;

    FunctionParser parser;
    std::string expression;
};

#endif /* FUNCTIONPARSER_HH_ */
