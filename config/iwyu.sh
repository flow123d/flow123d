#!/bin/bash

###########################################
# Installation of include-what-you-use
##########################################
# meant to be modified manually

set -x

# Assuming default clang version
sudo apt install clang libclang-dev

# GIT
SRCDIR=include-what-you-use
clang_ver=`clang --version |grep version| sed 's/clang version \([0-9].[0-9].[0-9]\).*/\1/'`
clang_ver_short=${clang_ver%.*}
cd $SRCDIR && git checkout "clang_${clang_ver_short}"
cd ..

# BUild dir
BUILD=IWYU_build
rm -rf $BUILD
mkdir $BUILD
cd $BUILD

export CC=clang
export CXX=clang++
cmake -G "Unix Makefiles" -DIWYU_LLVM_ROOT_PATH=/usr/lib/llvm-${clang_ver_short} ../$SRCDIR

make

# default Linking fails
# Need to add linking for stdc++ and m

clang  -fno-rtti -std=c++11   CMakeFiles/include-what-you-use.dir/iwyu.cc.o CMakeFiles/include-what-you-use.dir/iwyu_ast_util.cc.o CMakeFiles/include-what-you-use.dir/iwyu_cache.cc.o CMakeFiles/include-what-you-use.dir/iwyu_driver.cc.o CMakeFiles/include-what-you-use.dir/iwyu_getopt.cc.o CMakeFiles/include-what-you-use.dir/iwyu_globals.cc.o CMakeFiles/include-what-you-use.dir/iwyu_include_picker.cc.o CMakeFiles/include-what-you-use.dir/iwyu_lexer_utils.cc.o CMakeFiles/include-what-you-use.dir/iwyu_location_util.cc.o CMakeFiles/include-what-you-use.dir/iwyu_output.cc.o CMakeFiles/include-what-you-use.dir/iwyu_path_util.cc.o CMakeFiles/include-what-you-use.dir/iwyu_preprocessor.cc.o CMakeFiles/include-what-you-use.dir/iwyu_verrs.cc.o  -o include-what-you-use -rdynamic -lclangFrontend -lclangSerialization -lclangDriver -lclangParse -lclangSema -lclangAnalysis -lclangAST -lclangBasic -lclangEdit -lclangLex -lLLVMX86AsmParser -lLLVMX86CodeGen -lLLVMX86Desc -lLLVMX86AsmPrinter -lLLVMX86Info -lLLVMX86Utils -lLLVMCodeGen -lLLVMipo -lLLVMScalarOpts -lLLVMInstCombine -lLLVMTransformUtils -lLLVMTarget -lLLVMAnalysis -lLLVMOption -lLLVMMCDisassembler -lLLVMMCParser -lLLVMMC -lLLVMObject -lLLVMBitReader -lLLVMCore -lLLVMSupport -lpthread -lz -lcurses -lform -L/usr/lib -lstdc++ -ldl -lm

# Install
INSTAL_PREFIX=${HOME}/local/IWYU

make DESTDIR=${INSTAL_PREFIX}
CLANGBIN=`which clang`
CLANGINCLUDE=${CLANGBIN%bin/*}/lib/clang/${clang_ver}/include

IWYU_INCLUDE=${INSTAL_PREFIX}/usr/local/lib/${CLANGINCLUDE#*/lib/}
mkdir -p ${IWYU_INCLUDE}
cp ${CLANGINCLUDE}/* ${IWYU_INCLUDE}