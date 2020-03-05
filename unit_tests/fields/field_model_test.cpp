/*
 * field_model_test.cpp
 *
 *  Created on: Feb 26, 2020
 *      Author: David Flanderka
 *
 *  Tests EvalPoints, Integral classes ...
 */

#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <mesh_constructor.hh>
#include "arma_expect.hh"
#include <iostream>
#include <tuple>
#include <string>
#include <utility>
#include <type_traits>

#include "fields/field_model.hh"


typedef Field<3, FieldValue<3>::Scalar > ScalarField;
typedef Field<3, FieldValue<3>::VectorFixed > VectorField;


TEST(FieldModelTest, own_model) {
	ScalarField f_scal;
	VectorField f_vec;
	Fn functor;

	FieldModel<3, FieldValue<3>::VectorFixed, Fn, FieldValueCache<double>, FieldValueCache<double>> f_product(functor, f_scal.value_cache(), f_vec.value_cache());
	f_product.cache_update();
}

/*
 * Expand of std::tuple - functional solution:
 * https://stackoverflow.com/questions/10604794/convert-stdtuple-to-stdarray-c11
 */
// ------------- UTILITY---------------
/*#include <iostream>
#include <tuple>
#include <array>

template<int... Indices>
struct indices {
    using next = indices<Indices..., sizeof...(Indices)>;
};

template<int Size>
struct build_indices {
    using type = typename build_indices<Size - 1>::type::next;
};

template<>
struct build_indices<0> {
    using type = indices<>;
};

template<typename T>
using Bare = typename std::remove_cv<typename std::remove_reference<T>::type>::type;

template<typename Tuple>
constexpr
typename build_indices<std::tuple_size<Bare<Tuple>>::value>::type
make_indices()
{ return {}; }

template<typename Tuple, int... Indices>
std::array<
  typename std::tuple_element<0, Bare<Tuple>>::type,
    std::tuple_size<Bare<Tuple>>::value
>
to_array(Tuple&& tuple, indices<Indices...>)
{
    using std::get;
    return {{ get<Indices>(std::forward<Tuple>(tuple))... }};
}

template<typename Tuple>
auto to_array(Tuple&& tuple)
-> decltype( to_array(std::declval<Tuple>(), make_indices<Tuple>()) )
{
    return to_array(std::forward<Tuple>(tuple), make_indices<Tuple>());
}

TEST(FieldModelTest, tuple) {
    std::tuple<double, double, double, double> tup(1.5, 2.5, 4.5, 5.5);
    auto arr = to_array(tup);
    for (double x : arr)
        std::cout << x << " ";
    std::cout << std::endl;
} //*/

/*
 * Expand of std::tuple - needs -std=c++14:
 * https://genix.wordpress.com/2015/02/14/tuples-expanding-into-function-parameters/
 */
/*#if _MSC_VER == 1900 // hack for VS2015 CTP 5 broken decltype(auto) deduction in tuple_into_callable_n
#define TUPLE_FWD_RETURN(x) std::conditional_t< std::is_rvalue_reference<decltype(x)>::value, std::remove_reference_t<decltype(x)>, decltype(x)>(x)
#else
#define TUPLE_FWD_RETURN(x) x
#endif

// support for expanding tuples into an overloaded function call's arguments
#define wrap_overload(func) [](auto&&... ps){ return func( std::forward<decltype(ps)>(ps)... ); }

namespace detail
{
    //
    // base case for building up arguments for the function call
    //
    template< typename CALLABLE, typename TUPLE, int INDEX >
    struct tuple_into_callable_n
    {
        template< typename... Vs >
        static auto apply(CALLABLE f, TUPLE t, Vs&&... args) -> decltype(auto) // error: expected primary-expression before 'auto'
                                                                               // error: expected ')' before 'auto'
                                                                               // error: expected type-specifier before 'decltype'
                                                                               // error: expected initializer before 'decltype'
        {
            return tuple_into_callable_n<CALLABLE, TUPLE, INDEX - 1>::apply(
                f,
                std::forward<decltype(t)>(t),
                std::get<INDEX - 1>(std::forward<decltype(t)>(t)),
                std::forward<Vs>(args)...
            );
        }
    };

    //
    // terminal case - do the actual function call
    //
    template< typename CALLABLE, typename TUPLE >
    struct tuple_into_callable_n< CALLABLE, TUPLE, 0 >
    {
        template< typename... Vs >
        static auto apply(CALLABLE f, TUPLE t, Vs&&... args) -> decltype(auto) // same errors as previous
        {
            return TUPLE_FWD_RETURN(f(std::forward<Vs>(args)...));
        };
    };
}

template< typename FUNC, typename TUPLE >
auto tuple_into_callable(FUNC f, TUPLE&& t) -> decltype(auto) // same errors as previous
{
    return
        detail::tuple_into_callable_n<
            FUNC,
            decltype(t),
            std::tuple_size< std::remove_reference_t<TUPLE> >::value
        >::apply(f, std::forward<decltype(t)>(t) );
}

using namespace std;

class MoveOnly
{
public:
    MoveOnly(const MoveOnly&) = delete;
    MoveOnly(MoveOnly&& m) {};

    static MoveOnly Create()
    {
        return MoveOnly();
    }

private:
    MoveOnly() {}
};


string two_params(int i, int j)
{
    return "two_params int overload called.";
}

string two_params(float const& i, float const& j)
{
    return "two_params float const& overload called.";
}

string two_params(float&& i, float&& j)
{
    return "two_params float&& overload called.";
}

string move_only_receiver(MoveOnly&& m)
{
    return "move_only_receiver called.";
}

string no_params()
{
    return "no_params called.";
}

TEST(FieldModelTest, tuple) {

    // lvalue tuple
    auto t = make_tuple(1, 2);
    cout << tuple_into_callable(wrap_overload(two_params), t) << endl; // error: use of 'auto' in lambda parameter declaration only available with -std=c++14 or -std=gnu++14
                                                                       // error: expansion pattern 'int&&' contains no argument packs
                                                                       // error: 'ps' was not declared in this scope
                                                                       // error: 'tuple_into_callable' was not declared in this scope

    // rvalue tuple
    cout << tuple_into_callable(wrap_overload(two_params), make_tuple(1.0f, 2.0f)) << endl; // same errors as previous

    // const tuple
    auto const ct = make_tuple(1.0f, 2.0f);
    cout << tuple_into_callable(wrap_overload(two_params), ct) << endl; // same errors as previous

    // empty tuple -> empty function
    auto et = make_tuple();
    cout << tuple_into_callable(no_params, et) << endl;

    // tuple with move-only type
    auto move_only = MoveOnly::Create();
    auto mt = make_tuple(move(move_only));
    cout << tuple_into_callable(move_only_receiver, move(mt)) << endl; // note: tuple must be move'd into the callable so MoveOnly can be treated as &&
 }


/*
 * Expand of std::tuple - functional solution:
 * https://stackoverflow.com/questions/687490/how-do-i-expand-a-tuple-into-variadic-template-functions-arguments
 */
// ------------- UTILITY---------------
template<int...> struct index_tuple{};

template<int I, typename IndexTuple, typename... Types>
struct make_indexes_impl;

template<int I, int... Indexes, typename T, typename ... Types>
struct make_indexes_impl<I, index_tuple<Indexes...>, T, Types...>
{
    typedef typename make_indexes_impl<I + 1, index_tuple<Indexes..., I>, Types...>::type type;
};

template<int I, int... Indexes>
struct make_indexes_impl<I, index_tuple<Indexes...> >
{
    typedef index_tuple<Indexes...> type;
};

template<typename ... Types>
struct make_indexes : make_indexes_impl<0, index_tuple<>, Types...>
{};

// ----------UNPACK TUPLE AND APPLY TO FUNCTION ---------
template<class Ret, class... Args, int... Indexes >
Ret apply_helper( Ret (*pf)(Args...), index_tuple< Indexes... >, std::tuple<Args...>&& tup)
{
    return pf( std::forward<Args>( std::get<Indexes>(tup))... );
}

template<class Ret, class ... Args>
Ret apply(Ret (*pf)(Args...), const std::tuple<Args...>&  tup)
{
    return apply_helper(pf, typename make_indexes<Args...>::type(), std::tuple<Args...>(tup));
}

template<class Ret, class ... Args>
Ret apply(Ret (*pf)(Args...), std::tuple<Args...>&&  tup)
{
    return apply_helper(pf, typename make_indexes<Args...>::type(), std::forward<tuple<Args...>>(tup));
}

// --------------------- TEST ------------------
void one(int i, double d)
{
    std::cout << "function one(" << i << ", " << d << ");\n";
}
int two(int i)
{
    std::cout << "function two(" << i << ");\n";
    return i;
}
template<unsigned int size>
void vect(std::array<int, size> i_vec, std::array<double, size> d_vec)
{
    for (unsigned int i=0; i<size; ++i)
        std::cout << "Item " << i << ", function vect(" << i_vec[i] << ", " << d_vec[i] << ");\n";
}

TEST(FieldModelTest, tuple) {
    std::tuple<int, double> tup(23, 4.5);
    apply(one, tup);

    int d = apply(two, std::make_tuple(2));

    std::tuple<std::array<int, 3>, std::array<double, 3>> vec_tup(std::array<int, 3>({1, 2, 3}), std::array<double, 3>({4.5, 5.5, 6.5}) );
    apply(vect<3>, vec_tup);
} // */

// ----------------------------------------------------------------------------------------------------------------
// Pass of non-fixed size matrix as argument with fixed size.
template <uint nr, uint nc>
class MyMatrix {
public:
    typedef typename arma::mat::fixed<nr, nc> MatType;

    MyMatrix(double x) {
        uint n=1;
        for (uint i=0; i<nr; ++i)
            for (uint j=0; j<nc; ++j, ++n)
                matrix_(i,j) = n*x;
    }

    arma::mat matrix() const {
        return matrix_;
    }
private:
    MatType matrix_;
};

arma::vec3 product(typename arma::Col<double>::template fixed<1> a, arma::vec3 b) {
    return b * a;
}

TEST(FieldModelTest, armamats) {
    MyMatrix<1,1> scalar(0.5);
    MyMatrix<3,1> vector(2.0);

    auto res = product(scalar.matrix(), vector.matrix());
    std::cout << "Result = " << res;
}

