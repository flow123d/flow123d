/*
 * time_governor_test.cpp
 *
 *  Created on: May 20, 2011
 *      Author: jb
 */

#include <gtest/gtest.h>

#include "system/system.hh"
#include "coupling/time_marks.hh"
#include "coupling/time_governor.hh"

#include <input/input_type.hh>
#include <input/json_to_storage.hh>
#include <input/accessors.hh>

#include <limits>

#define EQUAL(a,b) INPUT_CHECK( (a) == (b), #a": %f and "#b":%f differs\n",a,b);


/**
 * In place input file.
 * CXX flags needs to be added for reading multiple line string
 * "-std=gnu++0x" (Use C++11 for auto and multiline raw strings)
 */
const string flow_json = R"JSON(
{
global_mesh = "some.msh",
equations= 
[
   {
      TYPE="EquationA",
      mesh={REF:"/global_mesh"},
      parameter_a=3.14 
   }  
],

time = { 
    start_time = 0.0, 
    end_time = 20.0
  }
}
)JSON";

/**
 * Test for class TimeMark only
 */
TEST (TimeMark, time_mark)
{
    TimeMark tm1(-1.0, 0x01);
    TimeMark tm2(0.0, TimeMark::every_type);
    TimeMark tm3(10.0,0x05);
    TimeMark::Type my_type = 0x0a;	
    
    //checking time values
    EQUAL(tm1.time(), -1.0);
    EQUAL(tm2.mark_type(), ~0x0);
    
    //checking type values - comparing masks
    EXPECT_TRUE(tm1.match_mask(0x01));
    EXPECT_TRUE(tm2.match_mask(~0x00));
    
    //adding type and checking mask
    tm3.add_to_type(my_type);
    EXPECT_TRUE(tm3.match_mask(0x0f));		//0x05 + 0x0a = 0x0f
    
    //comparing times
    EXPECT_TRUE( (tm1<tm2) && (tm2<tm3) );
}


/**
 * Test for class TimeMark, TimeMarks, TimeMarkIterator, TimeGovernor
 */
TEST (TimeGovernor, time_governor_marks_iterator)
{
    //creating shortcut to static time marks of time governor
    TimeMarks &tm = TimeGovernor::marks();
    
    //adding time marks
    tm.add(TimeMark(-1.0,tm.type_fixed_time() ) );
    tm.add(TimeMark(100,TimeMark::every_type));
    
    //creating mark types of our own
    xprintf(MsgDbg, "\nPredefined marktypes:\n\ttype_fixed_time_: %d\n\ttype_output_: %d\n\ttype_bc_change_: %d\n", tm.type_fixed_time(), tm.type_output(), tm.type_bc_change());
    // type_fixed_time_, type_output_ type_bc_change_ has been created -> 
    // new_mark_type = next_mark_type = 0x05<<1 = 0x08
    TimeMark::Type my_mark_type = tm.new_mark_type();	
    TimeMark::Type your_mark_type = 0x10;
    
    //adding equally spaced marks of my_mark_type type
    tm.add_time_marks(2.0, 0.25, 3.0, my_mark_type);
    tm.add(TimeMark(0.8, my_mark_type));
    tm.add(TimeMark(5.0, your_mark_type | tm.type_fixed_time() ));
    
    //inserting TimeMark
    tm.add(TimeMark(1.0,0x08));
    //adding mark type at previously defined time
    tm.add(TimeMark(3.0, your_mark_type));
    
    //constructing Time Governor from json input string
    
    static Input::Type::Record in_rec("RootInput", "Root record.");
    
    if (! in_rec.is_finished()) {
    in_rec.declare_key("time", TimeGovernor::input_type, "");
    in_rec.finish();
    }
    
    Input::JSONToStorage json_reader;
    //creating stream out of the json input string 
    std::stringstream in_stream(flow_json);
    //json reading according to keys defined in in_rec
    json_reader.read_stream(in_stream, in_rec);
    
    //getting root record
    Input::Record input = json_reader.get_root_interface<Input::Record>();
    
    //constructing time governor with time record read from input
    TimeGovernor *tm_tg = new TimeGovernor(  input.val<Input::Record>("time"), tm.type_fixed_time() );
    
    cout << tm;
    
    xprintf(MsgDbg, "\nTimeGovernor:\n\tstart_time: %f\n\tend_time: %f\n\tend_of_fixed_dt: %f\n\t time_step: %f\n\t last_time_step: %f\n", 
	    tm_tg->t(), tm_tg->end_time(),tm_tg->end_of_fixed_dt(), tm_tg->dt(), tm_tg->last_dt() );
	 
    //testing if TimeGovernor was correctly constructed
    EQUAL( tm_tg->t(), 0.0 );
    EQUAL( tm_tg->end_time(), 20.0 );
    EQUAL( tm_tg->end_of_fixed_dt(), 0.0);
    EQUAL( tm_tg->dt(), 20.0);
    EQUAL( tm_tg->last_dt(), 0.0);
    EQUAL( tm_tg->tlevel(), 0 );
    EXPECT_TRUE(tm_tg->is_changed_dt()); 	//changed from ZERO
	 
    //----------------- first testing of TimeMarks with TimeGovernor
    //next() is now done always from the time of governor = 0.0
    TimeMarks::iterator it = tm.begin();
    
    EQUAL( it->mark_type(), TimeMark::every_type ); //begins with ~0x0
    EQUAL( (++it)->time(), 100 );
    EQUAL( (--it)->time(), -numeric_limits<double>::infinity() );
    
    it = tm_tg->next(tm.type_fixed_time());
    EQUAL( (it)->time(), 5.0 );		//this goes from start time 0.0 not from -inf
    
    it = tm_tg->next(my_mark_type);
    EQUAL( (it)->time(), 0.8 );		//is included in every_type
    EQUAL( (++it)->time(), 1.0 );
    EQUAL( (++it)->time(), 2.0 );
    EQUAL( (++it)->time(), 2.25 );
    EQUAL( (++it)->time(), 2.5 );
    EQUAL( (++it)->time(), 2.75 );
    EQUAL( (++it)->time(), 3.0 );	//is included in 0x18
    
    it = tm_tg->next(your_mark_type);
    EQUAL( (it)->time(), 3.0 );	//is included in your_mark_type
    
    it = tm_tg->next(tm.type_fixed_time());
    EQUAL( (it)->time(), 5.0 );
    EQUAL( (++it)->time(), 20.0 );	//end_time = 20.0
    
    EQUAL( (--it)->time(), 5.0 );
    EQUAL( (--it)->time(), 0.0 );
    EQUAL( (--it)->time(), -1.0 );
    EQUAL( (--it)->time(), -numeric_limits<double>::infinity() );	//is included in every_type
    
    it = tm_tg->next( TimeMark::every_type);
    EQUAL( it->time(), 100 );
    EQUAL( (++it)->time(), numeric_limits<double>::infinity() );
    //-----------------
    
    //set permanent min_dt and max_dt
    tm_tg->set_permanent_constraint(0.01, 20.0);
    
    //testing setting of upper constraint
    //if out of allowed interval, cannot change the user constraints
    EQUAL( tm_tg->set_upper_constraint(25.0), -1);
    EQUAL( tm_tg->upper_constraint(), 20.0);
    EQUAL( tm_tg->set_upper_constraint(1e-4), 1);
    EQUAL( tm_tg->upper_constraint(), 20.0);
    
    //testing setting of lower constraint
    EQUAL( tm_tg->set_lower_constraint(25.0), -1);
    EQUAL( tm_tg->lower_constraint(), 0.01);
    EQUAL( tm_tg->set_lower_constraint(1e-4), 1);
    EQUAL( tm_tg->lower_constraint(), 0.01);
    
    //upper time step constraint fot next change of time_step
    EQUAL( tm_tg->set_upper_constraint(0.5), 0);
    
    //cout << "Estimated time (t+dt): " << tm_tg->estimate_time() << endl;
    //fixing time_step until next fixed mark_type (in time 5.0)
    tm_tg->fix_dt_until_mark();
    xprintf(MsgDbg, "Dt fixed. Estimated time (t+dt): %f", tm_tg->estimate_time() );
    
    
    //-----------------testing TimeGovernor's time stepping
    tm_tg->next_time();
    tm_tg->view();
   
    EQUAL( tm_tg->t(), 0.5 );
    EQUAL( tm_tg->tlevel(), 1 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_TRUE(tm_tg->is_changed_dt()); 	//changed from 20.0
    EXPECT_FALSE(tm_tg->is_current(TimeMark::every_type)); // time 0.5 (of tg) is not lt time 0.5 (of last timemark + dt = 0.0+0.5 = 0.5)
    //no time mark is in interval (0;0.5]
    
    tm_tg->next_time();
    tm_tg->view();

    EQUAL( tm_tg->t(), 1.0 );
    EQUAL( tm_tg->tlevel(), 2 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_FALSE(tm_tg->is_changed_dt()); 	//is still 0.5
    EXPECT_TRUE(tm_tg->is_current(my_mark_type)); // time 1.0 (of tg) is lt time 1.5 (of last timemark + dt = 1.0+0.5 = 1.5)
    // time mark 0x8 in time 1.0 is the last in interval (0.5;1.0]
    // time mark 0x8 in time 0.8 in interval (0.5;1.0] is skipped
    
    tm_tg->next_time();
    tm_tg->view();

    EQUAL( tm_tg->t(), 1.5 );
    EQUAL( tm_tg->tlevel(), 3 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_FALSE(tm_tg->is_changed_dt()); 	//is still 0.5
    EXPECT_FALSE(tm_tg->is_current(my_mark_type)); // time 1.5 (of tg) is NOT lt time 1.5 (of last timemark + dt = 1.0+0.5 = 1.5)
    //no time mark is in interval (1.0;1.5]
    
    tm_tg->next_time();
    tm_tg->view();

    EQUAL( tm_tg->t(), 2.0 );
    EQUAL( tm_tg->tlevel(), 4 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_FALSE(tm_tg->is_changed_dt()); 	//is still 0.5
    EXPECT_TRUE(tm_tg->is_current(my_mark_type)); // time 2.0 (of tg) is lt time 2.5 (of last timemark + dt = 2.0+0.5 = 2.5)
    // time mark 0x8 in time 2.0 is the last in interval (1.5;2.0]
    
    
    tm_tg->set_upper_constraint(2.0);
    //fixing time_step until next fixed mark_type (in time 5.0)
    tm_tg->fix_dt_until_mark();
    cout << "Dt fixed. Estimated time (t+dt): " << tm_tg->estimate_time() << endl;
    
    tm_tg->next_time();
    tm_tg->view();

    EQUAL( tm_tg->t(), 3.5 );
    EQUAL( tm_tg->tlevel(), 5 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_TRUE(tm_tg->is_changed_dt()); 	//is changed from 0.5 to 1.5
    EXPECT_FALSE(tm_tg->is_current(TimeMark::every_type)); // time 3.5 (of tg) is NOT lt time 1.5 (of last timemark + dt =0.0+1.5 = 1.5)
    EXPECT_TRUE(tm_tg->is_current(my_mark_type)); 
    EXPECT_TRUE(tm_tg->is_current(your_mark_type)); // time 3.5 (of tg) is lt time 4.5 (of last timemark + dt =3.0+1.5 = 4.5)
    // time mark 0x18 (0x08 included) in time 3.0 is the last in interval (2.0;3.5]
    // time mark 0x18 (0x10 included) in time 3.0 is the last in interval (2.0;3.5]
    
    tm_tg->next_time();
    tm_tg->view();

    EQUAL( tm_tg->t(), 5.0 );
    EQUAL( tm_tg->tlevel(), 6 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_FALSE(tm_tg->is_changed_dt()); 	//is still 1.5
    EXPECT_FALSE(tm_tg->is_current(my_mark_type)); // time 5.0 (of tg) is NOT lt time 4.5 (of last timemark + dt =3.0+1.5 = 4.5)
    //no time mark 0x08 is in interval (3.5;5.0]
    EXPECT_TRUE(tm_tg->is_current(your_mark_type)); // time 5.0 (of tg) is lt time 6.5 (of last timemark + dt =5.0+1.5 = 6.5)
    // time mark 0x11 (0x10 included) in time 3.0 is the last in interval (2.0;3.5]
     
    //-----------------
    
    
    //-----------------testing TimeMarks, TimeMarkIterator and last()
    it = tm_tg->next(TimeMark::every_type);
    EQUAL( (it)->time(), 100 );
    
    it = tm_tg->last(TimeMark::every_type);
    EQUAL( (it)->time(), -numeric_limits<double>::infinity() );
    
    it = tm_tg->last(my_mark_type);	//is included there
    EQUAL( (it)->time(), 3.0 );
    
    it = tm_tg->last(your_mark_type);	//is equal to time of TimeGovernor
    EQUAL( (it)->time(), 5.0 );
    
    it = tm_tg->last(tm.type_fixed_time());	//is included
    EQUAL( (it)->time(), 5.0 );
    
    //is_end()?
    
    delete tm_tg;
}


/**
 * Test for class steady TimeGovernor
 */
TEST (TimeGovernor, steady_time_governor)
{
    //DEFAULT CONSTRUCTOR
    TimeGovernor *steady_tg = new TimeGovernor();
    
    steady_tg->view("first_steady");
    
    EQUAL( steady_tg->t(), 0.0 );
    EQUAL( steady_tg->end_time(), numeric_limits<double>::infinity() );
    EQUAL( steady_tg->end_of_fixed_dt(), 0.0);
    EQUAL( steady_tg->dt(), numeric_limits<double>::infinity());
    EQUAL( steady_tg->last_dt(), 0.0);
    EQUAL( steady_tg->tlevel(), 0 );
    EXPECT_TRUE(steady_tg->is_changed_dt()); 	//changed from ZERO
    EQUAL( steady_tg->estimate_dt(), numeric_limits<double>::infinity());
    
    steady_tg->next_time();
    
    EQUAL( steady_tg->t(), numeric_limits<double>::infinity() );
    EQUAL( steady_tg->end_time(), numeric_limits<double>::infinity() );
    EQUAL( steady_tg->end_of_fixed_dt(), 0.0);
    EQUAL( steady_tg->dt(), 0.0);
    EQUAL( steady_tg->last_dt(), numeric_limits<double>::infinity());
    EQUAL( steady_tg->tlevel(), 1 );
    EXPECT_FALSE(steady_tg->is_changed_dt());  
    EQUAL( steady_tg->estimate_dt(), 0.0);
    
    EQUAL( steady_tg->fix_dt_until_mark(), 0.0);
    
    //CONSTRUCTOR with INITIAL TIME
    TimeGovernor *steady_tg_2 = new TimeGovernor(555.0);
    steady_tg_2->view("second_steady");
    
    //cout << steady_tg_2->marks();
    
    EQUAL( steady_tg_2->t(), 555.0 );
    EQUAL( steady_tg_2->end_time(), numeric_limits<double>::infinity() );
    EQUAL( steady_tg_2->end_of_fixed_dt(), 555.0);
    EQUAL( steady_tg_2->dt(), numeric_limits<double>::infinity());
    EQUAL( steady_tg_2->last_dt(), 0.0);
    EQUAL( steady_tg_2->tlevel(), 0 );
    EXPECT_TRUE(steady_tg_2->is_changed_dt()); 	//changed from ZERO
    EQUAL( steady_tg_2->estimate_dt(), numeric_limits<double>::infinity());
    
    steady_tg_2->next_time();
    
    EQUAL( steady_tg_2->t(), numeric_limits<double>::infinity() );
    EQUAL( steady_tg_2->end_time(), numeric_limits<double>::infinity() );
    EQUAL( steady_tg_2->end_of_fixed_dt(), 555.0);
    EQUAL( steady_tg_2->dt(), 0.0);
    EQUAL( steady_tg_2->last_dt(), numeric_limits<double>::infinity());
    EQUAL( steady_tg_2->tlevel(), 1 );
    EXPECT_FALSE(steady_tg_2->is_changed_dt()); 
    EQUAL( steady_tg_2->estimate_dt(), 0.0);
    
    delete steady_tg;
    delete steady_tg_2;
}
