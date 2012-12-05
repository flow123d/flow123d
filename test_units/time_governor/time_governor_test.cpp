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
   },
   {
      TYPE="EquationB",
      mesh={REF:"/global_mesh"},
      parameter_b=314,
      substances = [ "Rn", "Cs", "I", "C" ]
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
    TimeMarks tm;
    
    //adding time marks
    tm.add(TimeMark(-1.0,tm.type_fixed_time() ) );
    tm.add(TimeMark(0.0,TimeMark::every_type));
    
    //creating mark types of our own
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
    
    cout << tm;
    
    //constructing Time Governor from json input string
    
    static Input::Type::Record in_rec("RootInput", "Root record.");
    
    if (! in_rec.is_finished()) {
    in_rec.declare_key("time", TimeGovernor::get_input_type(), "");
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
    TimeGovernor *tm_tg = new TimeGovernor(  input.val<Input::Record>("time"), tm, tm.type_fixed_time() );

    cout << "TimeGovernor:\n\tstart_time: " << tm_tg->t()
	 << "\n\t end_time: " << tm_tg->end_time() 
	 << "\n\t end_of_fixed_dt: " << tm_tg->end_of_fixed_dt() 
	 << "\n\t time_step: " << tm_tg->dt() 
	 << "\n\t last_time_step: " << tm_tg->last_dt()<< endl;
	 
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
    TimeMarks::iterator it = tm_tg->marks().begin();	
    
    EQUAL( it->mark_type(), TimeMark::every_type ); //begins with ~0x0
    EQUAL( (++it)->time(), 0.0 );
    EQUAL( (--it)->time(), -numeric_limits<double>::infinity() );
    
    it = tm_tg->next(my_mark_type);
    EQUAL( (it)->time(), 0.8 );
    EQUAL( (++it)->time(), 1.0 );
    EQUAL( (++it)->time(), 2.0 );
    EQUAL( (++it)->time(), 2.25 );
    EQUAL( (++it)->time(), 2.5 );
    EQUAL( (++it)->time(), 2.75 );
    EQUAL( (++it)->time(), 3.0 );	//is included in 0x18
    
    it = tm_tg->next(your_mark_type);
    EQUAL( (it)->time(), 3.0 );
    
    it = tm_tg->next(tm.type_fixed_time());
    EQUAL( (it)->time(), 5.0 );
    EQUAL( (--it)->time(), 0.0 );	//is included in every_type
    EQUAL( (--it)->time(), -1.0 );
    EQUAL( (--it)->time(), -numeric_limits<double>::infinity() );	//is included in every_type
    
    it = tm_tg->next( TimeMark::every_type);
    EQUAL( it->time(), numeric_limits<double>::infinity() );
    //-----------------
    
    //set permanent min_dt and max_dt
    tm_tg->set_permanent_constrain(0.01, 20.0);
    //upper time_step constrain fot next change of time_step
    tm_tg->set_constrain(0.5);
    
    //cout << "Estimated time (t+dt): " << tm_tg->estimate_time() << endl;
    //fixing time_step until next fixed mark_type (in time 5.0)
    tm_tg->fix_dt_until_mark();
    //cout << "end_of_fixed_dt: " << tm_tg->end_of_fixed_dt() << endl;
    
    
    //-----------------testing TimeGovernor's time stepping
    tm_tg->next_time();
    cout << "TimeGovernor: time = " << tm_tg->t() 
    	   << "\ttime_step = "  << tm_tg->dt() << endl;
    EQUAL( tm_tg->t(), 0.5 );
    EQUAL( tm_tg->tlevel(), 1 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_TRUE(tm_tg->is_changed_dt()); 	//changed from 20.0
    EXPECT_FALSE(tm_tg->is_current(TimeMark::every_type)); // time 0.5 (of tg) is not lt time 0.5 (of last timemark + dt = 0.0+0.5 = 0.5)
    
    tm_tg->next_time();
    cout << "TimeGovernor: time = " << tm_tg->t() 
    	   << "\ttime_step = "  << tm_tg->dt() << endl;
    EQUAL( tm_tg->t(), 1.0 );
    EQUAL( tm_tg->tlevel(), 2 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_FALSE(tm_tg->is_changed_dt()); 	//is still 0.5
    EXPECT_TRUE(tm_tg->is_current(my_mark_type)); // time 1.0 (of tg) is lt time 1.3 (of last timemark + dt = 0.8+0.5 = 1.3)
    
    tm_tg->next_time();
    cout << "TimeGovernor: time = " << tm_tg->t() 
    	   << "\ttime_step = "  << tm_tg->dt() << endl;
    EQUAL( tm_tg->t(), 1.5 );
    EQUAL( tm_tg->tlevel(), 3 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_FALSE(tm_tg->is_changed_dt()); 	//is still 0.5
    EXPECT_FALSE(tm_tg->is_current(my_mark_type)); // time 1.5 (of tg) is NOT lt time 1.3 (of last timemark + dt = 1.0+0.5 = 1.5)
    
    
    tm_tg->next_time();
    cout << "TimeGovernor: time = " << tm_tg->t() 
    	   << "\ttime_step = "  << tm_tg->dt() << endl;
    EQUAL( tm_tg->t(), 2.0 );
    EQUAL( tm_tg->tlevel(), 4 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_FALSE(tm_tg->is_changed_dt()); 	//is still 0.5
    EXPECT_TRUE(tm_tg->is_current(my_mark_type)); // time 2.0 (of tg) is lt time 2.5 (of last timemark + dt = 2.0+0.5 = 2.5)
    
    
    tm_tg->set_constrain(2.0);
    //fixing time_step until next fixed mark_type (in time 5.0)
    tm_tg->fix_dt_until_mark();
    cout << "Dt fixed. Estimated time (t+dt): " << tm_tg->estimate_time() << endl;
    
    tm_tg->next_time();
    cout << "TimeGovernor: time = " << tm_tg->t() 
    	   << "\ttime_step = "  << tm_tg->dt() << endl;
    EQUAL( tm_tg->t(), 3.5 );
    EQUAL( tm_tg->tlevel(), 5 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_TRUE(tm_tg->is_changed_dt()); 	//is changed from 0.5 to 1.5
    EXPECT_TRUE(tm_tg->is_current(my_mark_type)); // time 3.5 (of tg) is lt time 3.5 (of last timemark + dt =3.0+1.5 = 4.5)
    
    tm_tg->next_time();
    cout << "TimeGovernor: time = " << tm_tg->t() 
    	   << "\ttime_step = "  << tm_tg->dt() << endl;
    EQUAL( tm_tg->t(), 5.0 );
    EQUAL( tm_tg->tlevel(), 6 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_FALSE(tm_tg->is_changed_dt()); 	//is still 1.5
    EXPECT_FALSE(tm_tg->is_current(my_mark_type)); // time 5.0 (of tg) is NOT lt time 4.5 (of last timemark + dt =3.0+1.5 = 4.5)

    //-----------------
    
    
    //-----------------testing TimeMarks, TimeMarkIterator and last()
    it = tm_tg->next(TimeMark::every_type);
    EQUAL( (it)->time(), numeric_limits<double>::infinity() );
    
    it = tm_tg->last(TimeMark::every_type);
    EQUAL( (it)->time(), 0.0 );
    
    it = tm_tg->last(my_mark_type);	//is included there
    EQUAL( (it)->time(), 3.0 );
    
    it = tm_tg->last(your_mark_type);	//is equal to time of TimeGovernor
    EQUAL( (it)->time(), 5.0 );
    
    it = tm_tg->last(tm.type_fixed_time());	//is included
    EQUAL( (it)->time(), 5.0 );
    
    delete tm_tg;
    
}

