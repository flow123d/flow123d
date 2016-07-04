/*
 * armadillo_setup.cc
 *
 *  Created on: Jul 4, 2016
 *      Author: jb
 */





#include <ostream>
#include <sstream>
#include "system/armadillo_setup.hh"
#include "system/exc_common.hh"
#include "system/logger.hh"
#include <armadillo>

using namespace std;

/**
 * Auxiliary exception in order to change standard frames_to_cut_.
 */
class ExcArmadillo : public ExceptionBase {
public:
    ExcArmadillo()
    {
        this->frames_to_cut_ = {"arma_stop"};
    }

    void print_info(std::ostringstream &out) const override {
        ::internal::ExcStream estream(out, *this);
        estream << "Armadillo Message: " << EI_Message::val << "\n";
    }
};

/**
 * Auxiliary output buffer to catch the error messages and report stacktrace.
 */
class ArmaStreamBuf : public std::stringbuf {
protected:
    /**
     * Override sync to throw on error message.
     */
    int sync() override;
};

int ArmaStreamBuf::sync() {
    if (this->str().find("error") !=  string::npos) {
        // Can not throw here, since any exception is cached somewhere between armadillo call point and this method.
        // Must write out the stack here.
        auto e = ExcArmadillo();
        e << EI_Message(this->str());
        _LOG( Logger::MsgType::message ) << e.what();

    }
    std::cout << this->str() << std::endl;
    this->str().clear();
    return 1;
}


void armadillo_setup()
{
    static ArmaStreamBuf stream_buf;
    static std::ostream arma_stream(&stream_buf);
    arma::set_stream_err1(arma_stream);
    arma::set_stream_err2(arma_stream);
}

