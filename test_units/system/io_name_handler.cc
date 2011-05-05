/*
 * io_name_handler.cc
 *
 *  Created on: May 4, 2011
 *      Author: honza
 */

#include"system/system.hh"
#include<iostream>
#include<petscerror.h>
#include<petsc.h>
#include<petscsys.h>
#include<gtest/gtest.h>

TEST(IONameHandler, InputDir) {
	std::string val = IONameHandler::get_instance()->get_input_file_name("${INPUT}/mesh.950.v01.msh");
	EXPECT_EQ("/home/honza/flow123dDev/1.6.5_build/test_units/input/mesh.950.v01.msh", val);
}
TEST(IONameHandler, OutputDir) {
	EXPECT_EQ("/output/mesh.950.v01.msh", IONameHandler::get_instance()->get_output_file_name("mesh.950.v01.msh"));
}
TEST(IONameHandler, AddPlaceholderItem) {
	std::string newString = "haha";
	IONameHandler::get_instance()->add_placeholder_item("${TEMP}", newString);
	EXPECT_EQ("/output/" + newString + "/nic.msh", IONameHandler::get_instance()->get_output_file_name("${TEMP}/nic.msh"));
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  int ierr;
  ierr = PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
  return RUN_ALL_TESTS();
}
