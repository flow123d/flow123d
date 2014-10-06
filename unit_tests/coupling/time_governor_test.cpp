/*
 * time_governor_test.cpp
 *
 *  Created on: May 20, 2011
 *      Author: jb
 */

#include <flow_gtest.hh>

#include "system/system.hh"
#include "coupling/time_marks.hh"
#include "coupling/time_governor.hh"

#include <input/input_type.hh>
#include <input/json_to_storage.hh>
#include <input/accessors.hh>

#include <limits>

//#define EXPECT_EQ(a,b) INPUT_CHECK( (a) == (b), #a": %f and "#b":%f differs\n",a,b);


/**
 * In place input file.
 * CXX flags needs to be added for reading multiple line string
 * "-std=gnu++0x" (Use C++11 for auto and multiline raw strings)
 */
const string flow_json = R"JSON(
{
time = { 
    start_time = 0.0, 
    end_time = 20.0
  }
}
)JSON";


const string json_with_init_dt = R"JSON(
{
time = { 
    start_time = 10.0, 
    end_time   = 30.0,
    init_dt    = 15
  }
}
)JSON";


const double inf_time = TimeGovernor::inf_time;



Input::Record read_input(const string &json_input)
{
	static Input::Type::Record in_rec("RootInput", "Root record.");

	if (! in_rec.is_finished()) {
	in_rec.declare_key("time", TimeGovernor::input_type, Input::Type::Default::obligatory(), "");
	in_rec.finish();
	}

	//json reading according to keys defined in in_rec
	Input::JSONToStorage json_reader(json_input, in_rec);

	//getting root record
	return json_reader.get_root_interface<Input::Record>().val<Input::Record>("time");
}



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
    EXPECT_EQ(tm1.time(), -1.0);
    EXPECT_EQ(tm2.mark_type(), TimeMark::Type(~0x0));
    
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
    // type_fixed_time_, type_output_ type_bc_change_ has been created -> 
    // new_mark_type = next_mark_type = 0x05<<1 = 0x08
    TimeMark::Type my_mark_type = tm.new_mark_type();	
    TimeMark::Type your_mark_type = tm.new_mark_type();
    
    //adding equally spaced marks of my_mark_type type
    tm.add_time_marks(2.0, 0.25, 3.0, my_mark_type);
    tm.add(TimeMark(0.8, my_mark_type));
    tm.add(TimeMark(5.0, your_mark_type | tm.type_fixed_time() ));
    
    //inserting TimeMark
    tm.add(TimeMark(1.0,my_mark_type | tm.type_input()  ));
    //adding mark type at previously defined time
    tm.add(TimeMark(3.0, your_mark_type));
    
    //constructing Time Governor from json input string
    
    static Input::Type::Record in_rec("RootInput", "Root record.");
    
    if (! in_rec.is_finished()) {
    in_rec.declare_key("time", TimeGovernor::input_type, Input::Type::Default::obligatory(), "");
    in_rec.finish();
    }
    
    //json reading according to keys defined in in_rec
    Input::JSONToStorage json_reader(flow_json, in_rec);
    
    //getting root record
    Input::Record input = json_reader.get_root_interface<Input::Record>();
    
    //constructing time governor with time record read from input
    TimeGovernor *tm_tg = new TimeGovernor(  input.val<Input::Record>("time"), my_mark_type  );
    
    cout << tm;
    
	 
    //testing if TimeGovernor was correctly constructed
    EXPECT_EQ( tm_tg->t(), 0.0 );
    EXPECT_EQ( tm_tg->end_time(), 20.0 );
    EXPECT_EQ( tm_tg->end_of_fixed_dt(), 0.0);
    EXPECT_EQ( tm_tg->dt(), 20.0);
    EXPECT_EQ( tm_tg->last_dt(), inf_time);
    EXPECT_EQ( tm_tg->tlevel(), 0 );
    EXPECT_TRUE(tm_tg->is_changed_dt()); 	//changed from ZERO
	 
    //----------------- first testing of TimeMarks with TimeGovernor
    //next() is now done always from the time of governor = 0.0
    TimeMarks::iterator it = tm.begin();
    
    EXPECT_EQ( it->mark_type(), TimeMark::every_type ); //begins with ~0x0
    EXPECT_EQ( (++it)->time(), 100 );
    EXPECT_EQ( (--it)->time(), -inf_time );
    
    it = tm_tg->next(tm.type_fixed_time());
    EXPECT_EQ( (it)->time(), 5.0 );		//this goes from start time 0.0 not from -inf
    
    it = tm_tg->next(my_mark_type);
    EXPECT_EQ( (it)->time(), 0.8 );		//is included in every_type
    EXPECT_EQ( (++it)->time(), 1.0 );
    EXPECT_EQ( (++it)->time(), 2.0 );
    EXPECT_EQ( (++it)->time(), 2.25 );
    EXPECT_EQ( (++it)->time(), 2.5 );
    EXPECT_EQ( (++it)->time(), 2.75 );
    EXPECT_EQ( (++it)->time(), 3.0 );	//is included in 0x18
    
    it = tm_tg->next(your_mark_type);
    EXPECT_EQ( (it)->time(), 3.0 );	//is included in your_mark_type
    
    it = tm_tg->next(tm.type_fixed_time());
    EXPECT_EQ( (it)->time(), 5.0 );
    EXPECT_EQ( (++it)->time(), 20.0 );	//end_time = 20.0
    
    EXPECT_EQ( (--it)->time(), 5.0 );
    EXPECT_EQ( (--it)->time(), 0.0 );
    EXPECT_EQ( (--it)->time(), -1.0 );
    EXPECT_EQ( (--it)->time(), -inf_time );	//is included in every_type
    
    it = tm_tg->next( TimeMark::every_type);
    EXPECT_EQ( it->time(), 100 );
    EXPECT_EQ( (++it)->time(), inf_time );
    //-----------------
    
    //set permanent min_dt and max_dt
    tm_tg->set_permanent_constraint(0.01, 20.0);
    
    //testing setting of upper constraint
    //if out of allowed interval, cannot change the user constraints
    EXPECT_EQ( tm_tg->set_upper_constraint(25.0), -1);
    EXPECT_EQ( tm_tg->upper_constraint(), 20.0);
    EXPECT_EQ( tm_tg->set_upper_constraint(1e-4), 1);
    EXPECT_EQ( tm_tg->upper_constraint(), 20.0);
    
    //testing setting of lower constraint
    EXPECT_EQ( tm_tg->set_lower_constraint(25.0), -1);
    EXPECT_EQ( tm_tg->lower_constraint(), 0.01);
    EXPECT_EQ( tm_tg->set_lower_constraint(1e-4), 1);
    EXPECT_EQ( tm_tg->lower_constraint(), 0.01);
    
    //upper time step constraint fot next change of time_step
    EXPECT_EQ( tm_tg->set_upper_constraint(0.5), 0);
    
    //cout << "Estimated time (t+dt): " << tm_tg->estimate_time() << endl;
    //fixing time_step until next fixed mark_type (in time 5.0)
    tm_tg->fix_dt_until_mark();
    //xprintf(MsgDbg, "Dt fixed. Estimated time (t+dt): %f", tm_tg->estimate_time() );
    
    
    //-----------------testing TimeGovernor's time stepping
    tm_tg->next_time();
    //tm_tg->view();
   
    EXPECT_EQ( tm_tg->t(), 0.5 );
    EXPECT_EQ( tm_tg->tlevel(), 1 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_TRUE(tm_tg->is_changed_dt()); 	//changed from 20.0
    EXPECT_FALSE(tm_tg->is_current(TimeMark::every_type)); // time 0.5 (of tg) is not lt time 0.5 (of last timemark + dt = 0.0+0.5 = 0.5)
    //no time mark is in interval (0;0.5]
    
    tm_tg->next_time();
    //tm_tg->view();

    EXPECT_EQ( tm_tg->t(), 1.0 );
    EXPECT_EQ( tm_tg->tlevel(), 2 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_FALSE(tm_tg->is_changed_dt()); 	//is still 0.5
    EXPECT_TRUE(tm_tg->is_current(my_mark_type) ); // time 1.0 (of tg) is lt time 1.5 (of last timemark + dt = 1.0+0.5 = 1.5)
    // time mark 0x8 in time 1.0 is the last in interval (0.5;1.0]
    // time mark 0x8 in time 0.8 in interval (0.5;1.0] is skipped
    
    tm_tg->next_time();
    //tm_tg->view();

    EXPECT_EQ( tm_tg->t(), 1.5 );
    EXPECT_EQ( tm_tg->tlevel(), 3 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_FALSE(tm_tg->is_changed_dt()); 	//is still 0.5
    EXPECT_FALSE(tm_tg->is_current(tm.type_fixed_time())); // time 1.5 (of tg) is NOT lt time 1.5 (of last timemark + dt = 1.0+0.5 = 1.5)
    //no time mark is in interval (1.0;1.5]
    
    tm_tg->next_time();
    //tm_tg->view();

    EXPECT_EQ( tm_tg->t(), 2.0 );
    EXPECT_EQ( tm_tg->tlevel(), 4 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_FALSE(tm_tg->is_changed_dt()); 	//is still 0.5
    EXPECT_TRUE(tm_tg->is_current(my_mark_type)); // time 2.0 (of tg) is lt time 2.5 (of last timemark + dt = 2.0+0.5 = 2.5)
    // time mark 0x8 in time 2.0 is the last in interval (1.5;2.0]
    
    
    tm_tg->set_upper_constraint(2.0);
    //fixing time_step until next fixed mark_type (in time 5.0)
    tm_tg->fix_dt_until_mark();
    cout << "Dt fixed. Estimated time (t+dt): " << tm_tg->estimate_time() << endl;
    
    tm_tg->next_time();
    //tm_tg->view();

    EXPECT_EQ( tm_tg->t(), 4.0 );
    EXPECT_EQ( tm_tg->tlevel(), 5 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_TRUE(tm_tg->is_changed_dt()); 	//is changed from 0.5 to 1.5
    EXPECT_FALSE(tm_tg->is_current(TimeMark::every_type)); // time 3.5 (of tg) is NOT lt time 1.5 (of last timemark + dt =0.0+1.5 = 1.5)
    EXPECT_TRUE(tm_tg->is_current(my_mark_type));
    EXPECT_TRUE(tm_tg->is_current(your_mark_type)); // time 3.5 (of tg) is lt time 4.5 (of last timemark + dt =3.0+1.5 = 4.5)
    // time mark 0x18 (0x08 included) in time 3.0 is the last in interval (2.0;3.5]
    // time mark 0x18 (0x10 included) in time 3.0 is the last in interval (2.0;3.5]
    
    tm_tg->next_time();
    //tm_tg->view();

    EXPECT_EQ( tm_tg->t(), 6.0 );
    EXPECT_EQ( tm_tg->tlevel(), 6 );
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_FALSE(tm_tg->is_changed_dt()); 	//is still 1.5
    EXPECT_FALSE(tm_tg->is_current(tm.type_fixed_time())); // time 5.0 (of tg) is NOT lt time 4.5 (of last timemark + dt =3.0+1.5 = 4.5)
    //no time mark 0x08 is in interval (3.5;5.0]
    EXPECT_FALSE(tm_tg->is_current(your_mark_type)); // time 5.0 (of tg) is lt time 6.5 (of last timemark + dt =5.0+1.5 = 6.5)
    // time mark 0x11 (0x10 included) in time 3.0 is the last in interval (2.0;3.5]
     
    //-----------------
    
    
    //-----------------testing TimeMarks, TimeMarkIterator and last()
    it = tm_tg->next(TimeMark::every_type);
    EXPECT_EQ( (it)->time(), 100 );
    
    it = tm_tg->last(TimeMark::every_type);
    EXPECT_EQ( (it)->time(), -inf_time );
    
    it = tm_tg->last(my_mark_type);	//is included there
    EXPECT_EQ( (it)->time(), 3.0 );
    
    it = tm_tg->last(your_mark_type);	//is equal to time of TimeGovernor
    EXPECT_EQ( (it)->time(), 5.0 );
    
    it = tm_tg->last(tm.type_fixed_time());	//is included
    EXPECT_EQ( (it)->time(), 5.0 );
    
    //is_end()?
    
    delete tm_tg;
}


TEST (TimeGovernor, simple_constructor)
{
	// Test of constructor without JSON input
	TimeGovernor tg(10, 0.5);

	EXPECT_EQ(tg.t(), 10);
	EXPECT_EQ(tg.dt(), 0.5);
	EXPECT_EQ(tg.end_time(), inf_time);
	EXPECT_EQ(tg.end_of_fixed_dt(), inf_time);

	tg.next_time();

	EXPECT_EQ(tg.t(), 10.5);

}


/**
 * Test for class steady TimeGovernor
 */
TEST (TimeGovernor, steady_time_governor)
{
    //DEFAULT CONSTRUCTOR
    TimeGovernor *steady_tg = new TimeGovernor();
    
    steady_tg->view("first_steady");
    
    EXPECT_EQ( steady_tg->t(), 0.0 );
    EXPECT_EQ( steady_tg->end_time(), inf_time );
    EXPECT_EQ( steady_tg->end_of_fixed_dt(), 0.0);
    EXPECT_EQ( steady_tg->dt(), inf_time);
    EXPECT_EQ( steady_tg->last_dt(), inf_time);
    EXPECT_EQ( steady_tg->tlevel(), 0 );
    EXPECT_TRUE(steady_tg->is_changed_dt()); 	//changed from ZERO
    EXPECT_EQ( steady_tg->estimate_dt(), 100);
    
    steady_tg->next_time();	// step to type_every mark at time 100
    
    EXPECT_EQ( steady_tg->t(), 100 );
    EXPECT_EQ( steady_tg->end_time(), inf_time );
    EXPECT_EQ( steady_tg->end_of_fixed_dt(), 0.0);
    EXPECT_EQ( steady_tg->dt(), 100);
    EXPECT_EQ( steady_tg->last_dt(), inf_time);
    EXPECT_EQ( steady_tg->tlevel(), 1 );
    EXPECT_TRUE(steady_tg->is_changed_dt());
    EXPECT_EQ( steady_tg->estimate_dt(), inf_time);
    
    EXPECT_EQ( steady_tg->fix_dt_until_mark(), 0.0);
    
    //CONSTRUCTOR with INITIAL TIME
    TimeGovernor *steady_tg_2 = new TimeGovernor(555.0);
    //steady_tg_2->view("second_steady");
    
    //cout << steady_tg_2->marks();
    
    EXPECT_EQ( steady_tg_2->t(), 555.0 );
    EXPECT_EQ( steady_tg_2->end_time(), inf_time );
    EXPECT_EQ( steady_tg_2->end_of_fixed_dt(), 555.0);
    EXPECT_EQ( steady_tg_2->dt(), inf_time);
    EXPECT_EQ( steady_tg_2->last_dt(), inf_time);
    EXPECT_EQ( steady_tg_2->tlevel(), 0 );
    EXPECT_TRUE(steady_tg_2->is_changed_dt()); 	//changed from ZERO
    EXPECT_EQ( steady_tg_2->estimate_dt(), inf_time);
    
    steady_tg_2->next_time();
    
    EXPECT_EQ( steady_tg_2->t(), inf_time );
    EXPECT_EQ( steady_tg_2->end_time(), inf_time );
    EXPECT_EQ( steady_tg_2->end_of_fixed_dt(), 555.0);
    EXPECT_EQ( steady_tg_2->dt(), inf_time);
    EXPECT_EQ( steady_tg_2->last_dt(), inf_time);
    EXPECT_EQ( steady_tg_2->tlevel(), 1 );
    EXPECT_FALSE(steady_tg_2->is_changed_dt()); 
    EXPECT_EQ( steady_tg_2->estimate_dt(), 0.0);
    
    delete steady_tg;
    delete steady_tg_2;
}
