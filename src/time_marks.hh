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
     *  through the TimeMarks::new_mark_type method.
     */
    typedef unsigned long int Type;

    /// Base mark type for strict time marks.
    static const Type strict;
    static const Type every_type;

    /**
     * Constructor for a TimeMarks::Mark
     * @param time - time of the mark
     * @param type - mark type
     *
     * <b> In order to create a strict TimeMark (at time=0.1) with base type output_type, use:
     * TimeMark( 0.1, output_type | TimeMark::strict)
     */
    TimeMark(double time, Type type) :
        time_(time), mark_type_(type) {}


    /// Getter for mark type.
    inline Type mark_type() const {
        return mark_type_;
    }

    /// True if the TimeMark is strict.
    inline bool is_strict() const {
        return mark_type_ & strict;
    }

    /// Getter for the time of the TimeMark.
    inline double time() const {
        return time_;
    }

    /**
     * Returns true if TimeMark's type has 1 on all positions where mask has 1.
     * @param mask {Select bits that should be 1 for matching mark types.
     */

    inline bool match_mask(const TimeMark::Type &mask) const {
        return ( mask & (~mark_type_) ) == 0;
    }

    /// Add more bits that a mark satisfies.
    inline void add_to_type(const TimeMark::Type &type) {
        mark_type_ |= type;
    }

    /// Comparison of time marks according to their time.
    bool operator<(const TimeMark& second) const
      { return time_ < second.time(); }


private:
    double time_;
    Type mark_type_;
};
/**
 * Output operator for TimeMark class.
 */
ostream& operator<<(ostream& stream, const TimeMark &marks);



/**
 * Iterator into the TimeMarks of particular mask. This is always const iterator, i.e. it points to const TimeMark.
 */
class TimeMarksIterator {
public:
    TimeMarksIterator(const vector<TimeMark> &marks,const  vector<TimeMark>::const_iterator &it, const TimeMark::Type &mask)
    : marks_(marks), it_(it), mask_(mask) {}

    TimeMarksIterator &operator=(const TimeMarksIterator &it)
    {ASSERT(&marks_ == &it.marks_, "Can not assign TimeMarks::iterator of different container.\n");
     it_=it.it_;
     mask_=it.mask_;
     return *this;
    }
    /// Prefix increment. Skip non-matching marks.
    TimeMarksIterator &operator++()
    { while ( it_ != marks_.begin() && ! (++it_) -> match_mask(mask_) ); return (*this); }

    /// Prefix decrement. Skip non-matching marks.
    TimeMarksIterator &operator--()
    { while ( it_ != marks_.end() && ! (--it_) -> match_mask(mask_) ); return (*this); }


    ///  * dereference operator
    inline const TimeMark & operator *() const
            { return *it_; }

    /// -> dereference operator
    inline const TimeMark * operator ->() const
            { return &(*(it_)); }

    inline bool operator ==(const TimeMarksIterator &other) const
        {return it_ == other.it_; }

    inline bool operator !=(const TimeMarksIterator &other) const
            {return it_ != other.it_; }

    TimeMark::Type mask()
    { return mask_; }
private:
    const vector<TimeMark> &marks_;
    vector<TimeMark>::const_iterator it_;
    TimeMark::Type mask_;
};

class TimeGovernor;

/**
 * @brief This class is a collection of time marks to manage various events occurring during simulation time.
 *
 * <b> TimeMark and their types </b>
 * One TimeMark consists of time and type (TimeMark::Type) see the constructor TimeMark::TimeMark.
 * The type of mark is bitmap where individual bits corresponds to some base event types like changing a BC, output solution, coupling time with another
 * equation and so on. Base types can be combined by bitwise or (operator|).
 *
 * There is one particular base mark type TimeMark::strict.
 * Only marks with this type are considered as fixed times by a TimeGovernor which is connected to particular TimeMarks object.
 *
 * <b> TimeMarks collection</b>
 * TimeMarks collect marks of various types and provides methods for iterating over stored marks. You can selectively access only marks matching given
 * type mask. See TimeMark::match_mask.
 *
 * You can add one new mark through method add or add evenly spaced marks of same type by TimeMarks::add_time_marks.
 *
 * You can allocate new TimeMark::Type in the context of one TimeMarks object by TimeMarks::new_mark_type and TimeMarks::new_strict_mark_type.
 *
 * For a given TimeGovernor (not necessarily connected one) you can ask about existence of mark in current time interval (TimeMarks::is_current) and iterate
 * around current time (TimeMarks::next and TimeMarks::last).
 *
 * In most cases there will be only one TimeMarks object for the whole solved problem and used by TimeGovernors of individual equations. However
 * this is not necessary.
 */
class TimeMarks {

public:
    /// Iterator class for iteration over time marks of particular type. This is always const_iterator.
    typedef TimeMarksIterator iterator;



    TimeMarks()
        : next_mark_type_(2)
    {
        marks_.push_back(TimeMark(-INFINITY, TimeMark::every_type));
        marks_.push_back(TimeMark(+INFINITY, TimeMark::every_type));
    }

    /**
     * Add a new base mark within the context of the particular TimeMarks instance.
     * User should keep the returned value (MarkType is basically a bitmap) for further queries and
     * TimeMark insertions. ATTENTION: You can not use the TimeMark::Type with other TimeMarks instance!
     */
    TimeMark::Type new_mark_type();
    /**
     * Same as @fn new_mark_type, but set the strict bit.
     * This distinguish those types that are used as a fix times for connected TimeGovernor.
     */
    TimeMark::Type new_strict_mark_type();

    /**
     * Basic method for inserting TimeMarks.
     * @param time    Time of the TimeMark.
     * @param type    MarkType or their combinations.
     */
    void add(const TimeMark &mark);

    /**
     * Method for creating and inserting equally spaced TimeMarks.
     * @param time    Time of the first TimeMark.
     * @param dt      Interval for further TimeMarks.
     * @param end_time  No marks after the end_time.
     * @param type    MarkType or their combinations.
     *
     * Current lazy implementation have complexity O(m*n) where m is number of inserted time merks and n number of time marks in the array.
     * TODO: O(n+m) implementation
     */
    void add_time_marks(double time, double dt, double end_time, TimeMark::Type type);

    /**
     * Find the last time mark matching given mask, and returns true if it is in the time interval of
     * current time step.
     */
    bool is_current(const TimeGovernor &tg, const TimeMark::Type &mask) const;

    /**
     * Return the first TimeMark with time strictly greater then tg.time() that match the mask.
     * The time governor tg  is used also for time comparisons.
     *
     * @param tg    the time governor
     * @param mask  mask of marks to iterate on
     *
     * TODO: have also method which accepts double (time) instead of the whole TimeGovernor.
     * and compare without safety.
     */
    TimeMarks::iterator next(const TimeGovernor &tg, const TimeMark::Type &mask) const;

    /**
     * Return the last TimeMark with time less or equal to tg.time() that match the mask.
     * The time governor tg  is used also for time comparisons.
     * @param tg    the time governor
     * @param mask  mask of marks to iterate on
     */
    TimeMarks::iterator last(const TimeGovernor &tg, const TimeMark::Type &mask) const;

    TimeMarks::iterator begin() const
    {return TimeMarks::iterator(marks_, marks_.begin(), TimeMark::every_type); }

    TimeMarks::iterator end() const
        {return TimeMarksIterator(marks_, --marks_.end(), TimeMark::every_type); }

    friend ostream& operator<<(ostream& stream, const TimeMarks &marks);

private:

    /// MarkType that will be used at next add_time_mark call.
    TimeMark::Type next_mark_type_;

    /// TimeMarks list sorted according to the their time.
    vector<TimeMark> marks_;
};

/**
 * Output operator for TimeMarks database.
 */
ostream& operator<<(ostream& stream, const  TimeMarks &marks);


#endif /* TIME_MARKS_HH_ */
