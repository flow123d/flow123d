***************************************************
** System of automatic tests in flow123d *********
***************************************************



********************
How to run tests?***
********************

To run tests, there is makefile system. Rules to run tests are: 

test - runs test specified by rule testbase in tests/makefile, default is 01_steady_flow_123d. To change, just modify rule testbase in tests/makefile

testall - run all tests specified by rule testall in tests/makefile

%.tst - run test named %, replace with test name. It's possible to wtite just start of name, for example - make 01.tst will run test 01_steady_flow_123d. If there is more tests found with same begining, returns Too many argument error.




***********************
How to add new test?***
***********************

To add new test, add dir with test into tests/ . You have to add makefile, if
you want to use makefile rules to run it, for example for tests
01_steady_flow_123d makefile should be like this :


INI_FILES="flow.ini"	// name of ini file
NPROC="1"	// number(s) of procs
FLOW_PARAMS=" -ksp_atol 1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor"	// aditional flow params

include ../all_tests.mk 


You should also modify tests/makefile and add your tests to rule makeall.




****************************
Where can I find outputs?***
****************************

All outputs are located in tests/Results/{inifile_name}.{nprocs} directory. All flow123d outputs, logs are here. And also log from comparison ref output and computed output from ndiff.pl script.




***********************************
How to use run_test.sh script?***
***********************************

You can also run tests without using makefile rules. Just run make_tests.sh script in tests/. Usage: 

./make_tests.sh {name_of_inifile} {nprocs} {flow_params}


Outputs are located in same dir as if you run it through makefiles.


If you have any other questions, email to: michal.nekvasil@tul.cz