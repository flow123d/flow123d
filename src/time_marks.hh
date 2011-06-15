/*
 * time_marks.hh
 *
 *  Created on: Jun 15, 2011
 *      Author: jb
 */

#ifndef TIME_MARKS_HH_
#define TIME_MARKS_HH_

/**
 * This class represents one record in the TimeMarks simple database.
 * Class members can not be modified after the item is created.
 */
class TimeMark {
public:
    /**
     *  MarkType is a bitmap where each bit represents one base type such as (strict, Output, Input, ...)
     *  This allow more complex queries through bitwise operations. Also one TimeMark can be shared by more events.
     *  In the context of TimeMarks the MarkType can be either strict or vague. If a TimeGovernor is connected to the TimeMarks object
     *  the  TimeMarks with strict MarkType are used to match exactly their times. Base MarkTypes should be obtained form TimeMarks class
     *  through the @fn new_mark_type method.
     */
    typedef unsigned long int Type;

    /// Base mark type for strict time marks.
    static const Type strict;
    /**
     * Constructor for a TimeMarks::Mark
     * @param time - time of the mark
     * @param type - mark type
     *
     * In order to create a strict TimeMark (at time=0.1) with base type output_type, use:
     * TimeMark( 0.1, output_type | TimeMark::strict)
     */
    TimeMark(double time, TimeMarkType type) :
        time_(time), mark_type_(type) {
    }

    /// Getter for mark type.
    TimeMarkType mark_type() const {
        return mark_type_;
    }

    /// True if the TimeMark is strict.
    bool is_strict() const {
        return mark_type_ && strict;
    }

    /// Getter for the time of the TimeMark.
    double time() const {
        return time_;
    }

    /// Comparison of time marks.
    bool operator<(const TimeMark &second_mark) {
        return time_ < second_mark.time();
    }

private:
    double time_;
    TimeMarkType mark_type_;
};
const TimeMark::Type TimeMark::strict = 0x1;

/**
 * Simple database of TimeMsrks. Provides questions about last and nearest TimeMarks for particular types. C
 */
class TimeMarks {

public:
    TimeMarks() :
        next_mark_type(2) {
    }

    /**
     * Add a new base mark within the context of the particular TimeMarks instance.
     * User should keep the returned value (MarkType is basically a bitmap) for further queries and
     * TimeMark insertions. ATTENTION: You can not use the TimeMark::Type with other TimeMarks instance!
     *
     * If @param strict is true, we set bit for strict types. This distinguish those types that are used as a fix times for connected TimeGovernor.
     */
    TimeMark::Type new_mark_type(bool strict) {
        ASSERT(next_mark_type != 0, "Can not allocate new mark type. The limit is 32 mark types.\n");
        TimeMark::Type current_type = next_mark_type;

        next_mark_type <<= 1;
        return current_type | TimeMark::strict;
    }

    /**
     * Basic method for inserting TimeMarks.
     * @par time    Time of the TimeMark.
     * @par type    MarkType or their combinations.
     */
    void add_time_mark(TimeMark mark) {
        vector<TimeMarks>::iterator first_ge = std::lower_bound(makrs.begin(), marks.end(), mark);
        marks.insert(first_ge, mark);
    }

    /**
     * Method for creating and inserting equally spaced TimeMarks.
     * @par time    Time of the first TimeMark.
     * @par dt      Interval for further TimeMarks.
     * @par end_time  No marks after the end_time.
     * @par type    MarkType or their combinations.
     *
     * Current lazy implementation have complexity O(m*n) where m is number of inserted time merks and n number of time marks in the array.
     * TODO: O(n+m) implementation
     */
    void add_time_marks(double time, double dt, double end_time, MarkType type) {
        for (double t = time; t < end_time; t += dt)
            add_time_mark(TimeMark(t, type));
    }

    bool is_current(const TimeGovernor &tg, const MarkType &mask) {
        TimeMark tm = last(tg, mask);
    }

    TimeMark &next(const TimeGovernor &tg, const MarkType &mask) {
        vector<TimeMarks>::iterator first_ge = std::lower_bound(makrs.begin(), marks.end(), mark);
        --first_ge;
        while (!tg.lt(first_ge.time()) || !mask_match(mask, first_ge.type()))
            ++first_ge;
        return first_ge;
    }

    /**
     * Return the last TimeMark before @param tg.time() that  match the @param mask.
     * The time governor @param tg  is used also for time comparisons.
     */
    TimeMark &last(const TimeGovernor &tg, const MarkType &mask) {
        vector<TimeMarks>::iterator first_ge = std::lower_bound(makrs.begin(), marks.end(), mark);
        while (!tg.ge(first_ge.time()) || !mask_match(mask, first_ge.type()))
            --first_ge;
        return first_ge;
    }

private:

    /// Returns true if @param type has 1 on all positions where @param mask has 1.
    bool mask_match(MarkType &mask, MarkType &type) {
        return mask & (~type) == 0;
    }

    /// MarkType that will be used at next add_time_mark call.
    MarkType next_mark_type;

    /// TimeMarks list sorted according to the their time.
    vector<TimeMark> marks;

    /// Since there can be queries form different TimeGovernors, we can not simply
    /// pass through marks. We rather iterate marks every time using the largest
    /// query time as the initial point.
    double guess_time;
    list<TimeMark> guess_iter;
};

#endif /* TIME_MARKS_HH_ */
