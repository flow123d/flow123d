/*
 * time_governor_test.cpp
 *
 *  Created on: May 20, 2011
 *      Author: jb
 */

#include <flow_gtest.hh>
#include "system/system.hh"
#include <input/input_type.hh>
#include <input/json_to_storage.hh>
#include <input/accessors.hh>
#include "tools/time_governor.hh"
#include "tools/time_marks.hh"






const double inf_time = TimeGovernor::inf_time;


/**
 * Auxiliary function to declare and read test input.
 */
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
	static Input::Record rec;
	rec = json_reader.get_root_interface<Input::Record>();
	return rec.val<Input::Record>("time");
}


TEST(TimeStep, all) {
    TimeStep step(2.0);
    EXPECT_EQ(2.0, step.end());
    EXPECT_EQ(1.0, step.length());
    EXPECT_EQ(0, step.index());

    TimeStep step1=step.make_next(1.5);
    EXPECT_EQ(3.5, step1.end());
    EXPECT_EQ(1.5, step1.length());
    EXPECT_EQ(1, step1.index());

    TimeStep step2=step1.make_next(1.0, 10 );
    EXPECT_EQ(10, step2.end());
    EXPECT_EQ(1.0, step2.length());
    EXPECT_EQ(2, step2.index());
}


TEST(TimeGovernor, step) {
    TimeGovernor::marks().reinit();
    string tg_in="{time = { start_time = 0.0, end_time = 10.0 } }";
    TimeGovernor tg( read_input(tg_in));
    tg.marks().add_time_marks(0.0, 1.0, 10.0, tg.equation_fixed_mark_type());
    EXPECT_EQ(0, tg.step().index());
    EXPECT_EQ(0, tg.step(-1).index());
    EXPECT_EQ(0, tg.step(0).index());

    EXPECT_THROW( {tg.step(1);}, TimeGovernor::ExcMissingTimeStep);
    EXPECT_THROW( {tg.step(-2);}, TimeGovernor::ExcMissingTimeStep);

    tg.next_time();
    EXPECT_EQ(1, tg.step().index());
    EXPECT_EQ(1, tg.step(-1).index());
    EXPECT_EQ(1, tg.step(1).index());
    EXPECT_EQ(0, tg.step(0).index());
    EXPECT_EQ(0, tg.step(-2).index());

    EXPECT_THROW( {tg.step(2);}, TimeGovernor::ExcMissingTimeStep);
    EXPECT_THROW( {tg.step(-3);}, TimeGovernor::ExcMissingTimeStep);

    tg.next_time();
    EXPECT_EQ(2, tg.step().index());
    EXPECT_EQ(2, tg.step(-1).index());
    EXPECT_EQ(2, tg.step(2).index());
    EXPECT_EQ(1, tg.step(1).index());
    EXPECT_EQ(1, tg.step(-2).index());

    EXPECT_THROW( {tg.step(0);}, TimeGovernor::ExcMissingTimeStep);
    EXPECT_THROW( {tg.step(3);}, TimeGovernor::ExcMissingTimeStep);

}

TEST(TimeGovernor, comparisons)
{
    TimeGovernor::marks().reinit();
    string tg_in="{time = { start_time = 0.0, end_time = 1E10 } }";
    TimeGovernor tg( read_input(tg_in));

    EXPECT_EQ(1.0, tg.dt());
    EXPECT_EQ(0.0, tg.t());

    // First we test around zero
    // test all comparators, later test just one, since other are correctly related
    EXPECT_TRUE( tg.le(0.0) );
    EXPECT_TRUE( tg.le(0.5) );
    EXPECT_TRUE( tg.le(-1E-15) );
    EXPECT_FALSE( tg.le(-0.1) );

    EXPECT_TRUE( tg.ge(0.0) );
    EXPECT_TRUE( tg.ge(-0.1) );
    EXPECT_TRUE( tg.ge(1E-15) );
    EXPECT_FALSE( tg.ge(0.1) );

    EXPECT_FALSE( tg.lt(0.0) );
    EXPECT_TRUE( tg.lt(0.5) );
    EXPECT_FALSE( tg.lt(1E-15) );
    EXPECT_FALSE( tg.lt(-0.1) );

    EXPECT_FALSE( tg.gt(0.0) );
    EXPECT_TRUE( tg.gt(-0.1) );
    EXPECT_FALSE( tg.gt(-1E-15) );
    EXPECT_FALSE( tg.gt(0.1) );

    double expect_dt=1.0;
    for(double exp_time=1.0;
            exp_time<tg.end_time();
            exp_time*=2) {
        tg.marks().add(TimeMark(exp_time, tg.equation_fixed_mark_type()));
        tg.next_time();
        //tg.view("eQ");


        EXPECT_EQ(expect_dt, tg.dt());
        expect_dt=exp_time;
        EXPECT_EQ(exp_time, tg.t());

        EXPECT_TRUE( tg.ge(exp_time) );
        EXPECT_TRUE( tg.ge(exp_time-1.0) );
        EXPECT_TRUE( tg.ge(exp_time*(1+1E-15)) );
        EXPECT_FALSE( tg.ge(exp_time+1.0) );
    }

}


TEST(TimeGovernor, estimate_dt)
{
    TimeGovernor::marks().reinit();
    string tg_in="{time = { start_time = 0.0, end_time = 100.0 } }";
    TimeGovernor tg( read_input(tg_in));
    tg.marks().add(TimeMark(0.5, tg.equation_fixed_mark_type()));
    EXPECT_FLOAT_EQ(0.5, tg.estimate_dt());
    tg.set_upper_constraint(0.15);
    EXPECT_FLOAT_EQ(0.125, tg.estimate_dt());
    tg.next_time();
    EXPECT_FLOAT_EQ(0.375, tg.estimate_dt());
    tg.next_time();

    // test slight rounding up, violating upper constraint
    tg.marks().add(TimeMark(1.0, tg.equation_fixed_mark_type()));
    tg.set_upper_constraint(0.1+1E-14);
    EXPECT_FLOAT_EQ(0.1, tg.estimate_dt());
    double epsilon = 0.4*numeric_limits<double>::epsilon();
    EXPECT_FALSE( 0.1 == (0.1-epsilon) );
    tg.set_upper_constraint(0.1-epsilon);
    EXPECT_EQ(0.1, tg.estimate_dt());
    tg.set_upper_constraint(0.1 - 2*epsilon);
    EXPECT_TRUE(0.1 > tg.estimate_dt());
}




/**
 * Test for class TimeMark, TimeMarks, TimeMarkIterator, TimeGovernor
 */
TEST (TimeGovernor, time_governor_marks_iterator)
{
    TimeGovernor::marks().reinit();
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
    
    cout << tm;
    const string flow_json = R"JSON(
    {
    time = { 
        start_time = 0.0, 
        end_time = 20.0
      }
    }
    )JSON";
    //constructing Time Governor from json input string
    TimeGovernor *tm_tg = new TimeGovernor( read_input(flow_json), my_mark_type  );
    
    cout << tm;
    
	 
    //testing if TimeGovernor was correctly constructed
    EXPECT_EQ(0.0 ,tm_tg->t());
    EXPECT_EQ(20.0 ,tm_tg->end_time());
    EXPECT_FLOAT_EQ( 0.0, tm_tg->end_of_fixed_dt());
//    EXPECT_EQ(20.0,tm_tg->dt());                                  // FAILURE
    EXPECT_EQ(inf_time,tm_tg->last_dt());
    EXPECT_EQ(0 ,tm_tg->tlevel());
    EXPECT_TRUE(tm_tg->is_changed_dt()); 	//changed from ZERO
	 
    //----------------- first testing of TimeMarks with TimeGovernor
    //next() is now done always from the time of governor = 0.0
    TimeMarks::iterator it = tm.begin();
    

    EXPECT_EQ(TimeMark::every_type ,it->mark_type()); //begins with ~0x0
    EXPECT_EQ(100 ,(++it)->time());
    EXPECT_EQ(-inf_time ,(--it)->time());
    
    it = tm_tg->next(tm.type_fixed_time());
    EXPECT_EQ(5.0 ,(it)->time());		//this goes from start time 0.0 not from -inf
    
    it = tm_tg->next(my_mark_type);
    EXPECT_EQ(0.8 ,(it)->time());		//is included in every_type
    EXPECT_EQ(1.0 ,(++it)->time());
    EXPECT_EQ(2.0 ,(++it)->time());
    EXPECT_EQ(2.25 ,(++it)->time());
    EXPECT_EQ(2.5 ,(++it)->time());
    EXPECT_EQ(2.75 ,(++it)->time());
    EXPECT_EQ(3.0 ,(++it)->time());	//is included in 0x18
    
    it = tm_tg->next(your_mark_type);
    EXPECT_EQ(3.0 ,(it)->time());	//is included in your_mark_type
    
    it = tm_tg->next(tm.type_fixed_time());
    EXPECT_EQ(5.0 ,(it)->time());
    EXPECT_EQ(20.0 ,(++it)->time());	//end_time = 20.0
    
    EXPECT_EQ(5.0 ,(--it)->time());
    EXPECT_FLOAT_EQ( 0.0, (--it)->time());
    EXPECT_EQ(-1.0 ,(--it)->time());
    EXPECT_EQ(-inf_time ,(--it)->time());	//is included in every_type
    
    it = tm_tg->next( TimeMark::every_type);
    EXPECT_EQ(100 ,it->time());
    EXPECT_EQ(inf_time ,(++it)->time());
    //-----------------
    
    //set permanent min_dt and max_dt
    tm_tg->set_permanent_constraint(0.01, 20.0);
    
    //testing setting of upper constraint
    //if out of allowed interval, cannot change the user constraints
    EXPECT_EQ(-1,tm_tg->set_upper_constraint(25.0));
    EXPECT_EQ(20.0,tm_tg->upper_constraint());
    EXPECT_EQ(1,tm_tg->set_upper_constraint(1e-4));
    EXPECT_EQ(20.0,tm_tg->upper_constraint());
    
    //testing setting of lower constraint
    EXPECT_EQ(-1,tm_tg->set_lower_constraint(25.0));
    EXPECT_EQ(0.01,tm_tg->lower_constraint());
    EXPECT_EQ(1,tm_tg->set_lower_constraint(1e-4));
    EXPECT_EQ(0.01,tm_tg->lower_constraint());
    
    //upper time step constraint fot next change of time_step
    EXPECT_EQ(0,tm_tg->set_upper_constraint(0.5));
    
    //cout << "Estimated time (t+dt): " << tm_tg->estimate_time() << endl;
    //fixing time_step until next fixed mark_type (in time 5.0)
    tm_tg->fix_dt_until_mark();
    //xprintf(MsgDbg, "Dt fixed. Estimated time (t+dt): %f", tm_tg->estimate_time() );
    
    //-----------------testing TimeGovernor's time stepping
    tm_tg->next_time();
    //tm_tg->view();
   
    EXPECT_EQ(0.5 ,tm_tg->dt());
    EXPECT_EQ(0.5 ,tm_tg->t());
    EXPECT_EQ(1 ,tm_tg->tlevel());
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_TRUE(tm_tg->is_changed_dt()); 	//changed to 20.0     // FAILURE

    // FAILURE
    EXPECT_FALSE(tm_tg->is_current(TimeMark::every_type)); // time 0.5 (of tg) is not lt time 0.5 (of last timemark + dt = 0.0+0.5 = 0.5)
    //no time mark is in interval (0;0.5]
    
    tm_tg->next_time();
    //tm_tg->view();

    EXPECT_EQ(1.0 ,tm_tg->t());
    EXPECT_EQ(2 ,tm_tg->tlevel());
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_FALSE(tm_tg->is_changed_dt()); 	//is still 0.5
    EXPECT_TRUE(tm_tg->is_current(my_mark_type) ); // time 1.0 (of tg) is lt time 1.5 (of last timemark + dt = 1.0+0.5 = 1.5)
    // time mark 0x8 in time 1.0 is the last in interval (0.5;1.0]
    // time mark 0x8 in time 0.8 in interval (0.5;1.0] is skipped
    
    tm_tg->next_time();
    //tm_tg->view();

    EXPECT_EQ(1.5 ,tm_tg->t());
    EXPECT_EQ(3 ,tm_tg->tlevel());
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_FALSE(tm_tg->is_changed_dt()); 	//is still 0.5

    // FAILURE
    EXPECT_FALSE(tm_tg->is_current(tm.type_fixed_time())); // time 1.5 (of tg) is NOT lt time 1.5 (of last timemark + dt = 1.0+0.5 = 1.5)
    //no time mark is in interval (1.0;1.5]
    
    tm_tg->next_time();
    //tm_tg->view();

    EXPECT_EQ(2.0 ,tm_tg->t());
    EXPECT_EQ(4 ,tm_tg->tlevel());
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

    EXPECT_EQ(4.0 ,tm_tg->t());
    EXPECT_EQ(5 ,tm_tg->tlevel());
    EXPECT_FALSE(tm_tg->is_end());
    // FAILURE
    EXPECT_TRUE(tm_tg->is_changed_dt()); 	//is changed from 0.5 to 1.5
    // FAILURE
    EXPECT_FALSE(tm_tg->is_current(TimeMark::every_type)); // time 3.5 (of tg) is NOT lt time 1.5 (of last timemark + dt =0.0+1.5 = 1.5)
    // FAILURE
    EXPECT_TRUE(tm_tg->is_current(my_mark_type));
    // FAILURE
    EXPECT_TRUE(tm_tg->is_current(your_mark_type)); // time 3.5 (of tg) is lt time 4.5 (of last timemark + dt =3.0+1.5 = 4.5)
    // time mark 0x18 (0x08 included) in time 3.0 is the last in interval (2.0;3.5]
    // time mark 0x18 (0x10 included) in time 3.0 is the last in interval (2.0;3.5]
    
    tm_tg->next_time();
    //tm_tg->view();

    EXPECT_EQ(6.0 ,tm_tg->t());
    EXPECT_EQ(6 ,tm_tg->tlevel());
    EXPECT_FALSE(tm_tg->is_end());
    EXPECT_FALSE(tm_tg->is_changed_dt()); 	//is still 1.5
    // FAILURE
    EXPECT_FALSE(tm_tg->is_current(tm.type_fixed_time())); // time 5.0 (of tg) is NOT lt time 4.5 (of last timemark + dt =3.0+1.5 = 4.5)
    //no time mark 0x08 is in interval (3.5;5.0]
    // FAILURE
    EXPECT_FALSE(tm_tg->is_current(your_mark_type)); // time 5.0 (of tg) is lt time 6.5 (of last timemark + dt =5.0+1.5 = 6.5)
    // time mark 0x11 (0x10 included) in time 3.0 is the last in interval (2.0;3.5]
     
    //-----------------
    
    
    //-----------------testing TimeMarks, TimeMarkIterator and last()
    it = tm_tg->next(TimeMark::every_type);
    EXPECT_EQ(100 ,(it)->time());
    
    it = tm_tg->last(TimeMark::every_type);
    EXPECT_EQ(-inf_time ,(it)->time());
    
    it = tm_tg->last(my_mark_type);	//is included there
    EXPECT_EQ(3.0 ,(it)->time());
    
    it = tm_tg->last(your_mark_type);	//is equal to time of TimeGovernor
    EXPECT_EQ(5.0 ,(it)->time());
    
    it = tm_tg->last(tm.type_fixed_time());	//is included
    EXPECT_EQ(5.0 ,(it)->time());
    
    delete tm_tg;
}


TEST (TimeGovernor, simple_constructor)
{
    TimeGovernor::marks().reinit();

	// Test of constructor without JSON input
	TimeGovernor tg(10, 0.5);

	EXPECT_EQ(tg.t(), 10);
//	EXPECT_EQ(tg.dt(), 0.5);
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
    TimeGovernor::marks().reinit();

    //DEFAULT CONSTRUCTOR
    TimeGovernor *steady_tg = new TimeGovernor();
    TimeMarks &tm=steady_tg->marks();
    tm.add(TimeMark(100,TimeMark::every_type));

    
    steady_tg->view("first_steady");
    tm.add(TimeMark(0.0, tm.type_output() | steady_tg->equation_mark_type()) );
    
    EXPECT_EQ( steady_tg->t(), 0.0 );
    EXPECT_EQ( steady_tg->end_time(), inf_time );
    EXPECT_FLOAT_EQ( 0.0, steady_tg->end_of_fixed_dt() );
//    EXPECT_EQ( inf_time, steady_tg->dt());
    EXPECT_EQ( steady_tg->last_dt(), inf_time);
    EXPECT_EQ( steady_tg->tlevel(), 0 );
    EXPECT_TRUE(steady_tg->is_current(TimeGovernor::marks().type_output()));

    EXPECT_TRUE(steady_tg->is_changed_dt()); 	//changed from ZERO
    EXPECT_EQ( steady_tg->estimate_dt(), 100);
    
    steady_tg->next_time();	// step to type_every mark at time 100
    
    EXPECT_EQ( steady_tg->t(), 100 );
    EXPECT_EQ( steady_tg->end_time(), inf_time );
    EXPECT_EQ( steady_tg->end_of_fixed_dt(), 0.0);
    // FAILURE
    EXPECT_EQ( 100, steady_tg->dt());
    EXPECT_EQ( steady_tg->last_dt(), 1.0);
    EXPECT_EQ( steady_tg->tlevel(), 1 );
    // FAILURE
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
    //EXPECT_EQ( steady_tg_2->dt(), inf_time);
    EXPECT_EQ( steady_tg_2->last_dt(), inf_time);
    EXPECT_EQ( steady_tg_2->tlevel(), 0 );
    EXPECT_TRUE(steady_tg_2->is_changed_dt()); 	//changed from ZERO
    EXPECT_EQ( steady_tg_2->estimate_dt(), inf_time);
    
    steady_tg_2->next_time();
    
    EXPECT_EQ( steady_tg_2->t(), inf_time );
    EXPECT_EQ( steady_tg_2->end_time(), inf_time );
    EXPECT_EQ( steady_tg_2->end_of_fixed_dt(), 555.0);
    EXPECT_EQ( steady_tg_2->dt(), inf_time);
    EXPECT_EQ( steady_tg_2->last_dt(), 1.0);
    EXPECT_EQ( steady_tg_2->tlevel(), 1 );
    EXPECT_TRUE(steady_tg_2->is_changed_dt());
    EXPECT_EQ( steady_tg_2->estimate_dt(), 0.0);
    
    delete steady_tg;
    delete steady_tg_2;
}
