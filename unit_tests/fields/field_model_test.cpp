/*
 * field_model_test.cpp
 *
 *  Created on: Feb 26, 2020
 *      Author: David Flanderka
 *
 *  Tests FieldModel class and expand tuple to arguments.
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
#include "mesh/accessors.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"


typedef Field<3, FieldValue<3>::Scalar > ScalarField;
typedef Field<3, FieldValue<3>::VectorFixed > VectorField;


using Scalar = typename arma::Col<double>::template fixed<1>;
using Vector = arma::vec3;
using Tensor = arma::mat33;


class FnProduct {
public:
    typedef std::tuple<
            Field<3,FieldValue<3>::Scalar>,
            Field<3, FieldValue<3>::VectorFixed>
    > DepFields;

    Vector operator() (Scalar a, Vector v) {
        return a(0) * v;
    }
};



class FnOther {
public:
    typedef std::tuple<
            Field<3, FieldValue<3>::VectorFixed>,
            Field<3,FieldValue<3>::Scalar>,
            Field<3, FieldValue<3>::VectorFixed>
    > DepFields;

    Vector operator() (Vector a, Scalar c, Vector b) {
        return a + c(0) * b;
    }
};


// Test of FieldModel - used objects and functionalities in field_model.hh.
TEST(FieldModelTest, own_model) {
    unsigned int n_items = 10; // number of tested items

    ScalarField f_scal;
    VectorField f_vec;
    FieldValueCache<typename FieldValue<3>::VectorFixed::element_type> fvc(FieldValue<3>::VectorFixed::NRows_, FieldValue<3>::VectorFixed::NCols_);
    std::vector< ElementAccessor<3> > element_set;

    // initialize field caches
    std::shared_ptr<EvalPoints> eval_points = std::make_shared<EvalPoints>();
    Quadrature *q_bulk = new QGauss(3, 2);
    Quadrature *q_side = new QGauss(2, 2);
    eval_points->add_bulk<3>(*q_bulk );
    eval_points->add_edge<3>(*q_side );
    fvc.init(eval_points, ElementCacheMap::n_cached_elements);
    f_scal.cache_allocate(eval_points);
    f_vec.cache_allocate(eval_points);

    // fill field caches
    arma::mat::fixed<1,1> scalar_val;
    arma::mat::fixed<3,1> vector_val;
    for (unsigned int i=0; i<n_items; ++i) {
        scalar_val(0,0) = 1.0 + i*0.5;
        vector_val(0,0) = 1.5 + 2*i;
        vector_val(1,0) = i + 0.1;
        vector_val(2,0) = 0.5 + i%2;
        f_scal.value_cache().data().set(i) = scalar_val;
        f_vec.value_cache().data().set(i) = vector_val;
    }

    {
        // FieldModel scalar * vector
        std::vector<arma::vec3> expected_vals = {{  1.50,  0.10, 0.50},
                                                 {  5.25,  1.65, 2.25},
                                                 { 11.00,  4.20, 1.00},
                                                 { 18.75,  7.75, 3.75},
                                                 { 28.50, 12.30, 1.50},
                                                 { 40.25, 17.85, 5.25},
                                                 { 54.00, 24.40, 2.00},
                                                 { 69.75, 31.95, 6.75},
                                                 { 87.50, 40.50, 2.50},
                                                 {107.25, 50.05, 8.25}};

    	auto f_product = Model<3, FieldValue<3>::VectorFixed>::create(FnProduct(), f_scal, f_vec);
        f_product.cache_update(fvc, 0, fvc.size(), element_set);
        for (unsigned int i=0; i<n_items; ++i) {
            auto val = fvc.data().template mat<3, 1>(i);
            EXPECT_ARMA_EQ(val, expected_vals[i]);
        }
    }

    {
        // FieldModel vector + scalar * vector
        std::vector<arma::vec3> expected_vals = {{  3.00,  0.20, 1.00},
                                                 {  8.75,  2.75, 3.75},
                                                 { 16.50,  6.30, 1.50},
                                                 { 26.25, 10.85, 5.25},
                                                 { 38.00, 16.40, 2.00},
                                                 { 51.75, 22.95, 6.75},
                                                 { 67.50, 30.50, 2.50},
                                                 { 85.25, 39.05, 8.25},
                                                 {105.00, 48.60, 3.00},
                                                 {126.75, 59.15, 9.75}};

        auto f_other = Model<3, FieldValue<3>::VectorFixed>::create(FnOther(), f_vec, f_scal, f_vec);
        f_other.cache_update(fvc, 0, fvc.size(), element_set);
        for (unsigned int i=0; i<n_items; ++i) {
            auto val = fvc.data().template mat<3, 1>(i);
            EXPECT_ARMA_EQ(val, expected_vals[i]);
        }
    }

}

// Following blocks are different auxiliary solutions only of development, partly will be modified and partly will be removed in final version of test.
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

//TEST(FieldModelTest, tuple) {
//
//    // lvalue tuple
//    auto t = make_tuple(1, 2);
//    cout << tuple_into_callable(wrap_overload(two_params), t) << endl;
//
//    // rvalue tuple
//    cout << tuple_into_callable(wrap_overload(two_params), make_tuple(1.0f, 2.0f)) << endl;
//
//    // const tuple
//    auto const ct = make_tuple(1.0f, 2.0f);
//    cout << tuple_into_callable(wrap_overload(two_params), ct) << endl;
//
//    // empty tuple -> empty function
//    auto et = make_tuple();
//    cout << tuple_into_callable(no_params, et) << endl;
//
//    // tuple with move-only type
//    auto move_only = MoveOnly::Create();
//    auto mt = make_tuple(move(move_only));
//    cout << tuple_into_callable(move_only_receiver, move(mt)) << endl; // note: tuple must be move'd into the callable so MoveOnly can be treated as &&
//} // */


/*
 * Expand of std::tuple - functional solution:
 * https://stackoverflow.com/questions/687490/how-do-i-expand-a-tuple-into-variadic-template-functions-arguments
 */
// ------------- UTILITY---------------
/*template<int...> struct index_tuple_test{};

template<int I, typename IndexTuple, typename... Types>
struct make_indexes_test_impl;

template<int I, int... Indexes, typename T, typename ... Types>
struct make_indexes_test_impl<I, index_tuple_test<Indexes...>, T, Types...>
{
    typedef typename make_indexes_test_impl<I + 1, index_tuple_test<Indexes..., I>, Types...>::type type;
};

template<int I, int... Indexes>
struct make_indexes_test_impl<I, index_tuple_test<Indexes...> >
{
    typedef index_tuple_test<Indexes...> type;
};

template<typename ... Types>
struct make_indexes_test : make_indexes_test_impl<0, index_tuple_test<>, Types...>
{};

// ----------UNPACK TUPLE AND APPLY TO FUNCTION ---------
template<class Ret, class... Args, int... Indexes >
Ret apply_helper_test( Ret (*pf)(Args...), index_tuple_test< Indexes... >, std::tuple<Args...>&& tup)
{
    return pf( std::forward<Args>( std::get<Indexes>(tup))... );
}

template<class Ret, class ... Args>
Ret apply_test(Ret (*pf)(Args...), const std::tuple<Args...>&  tup)
{
    return apply_helper_test(pf, typename make_indexes_test<Args...>::type(), std::tuple<Args...>(tup));
}

template<class Ret, class ... Args>
Ret apply_test(Ret (*pf)(Args...), std::tuple<Args...>&&  tup)
{
    return apply_helper_test(pf, typename make_indexes_test<Args...>::type(), std::forward<tuple<Args...>>(tup));
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
    apply_test(one, tup);

    int d = apply_test(two, std::make_tuple(2));

    std::tuple<std::array<int, 3>, std::array<double, 3>> vec_tup(std::array<int, 3>({1, 2, 3}), std::array<double, 3>({4.5, 5.5, 6.5}) );
    apply_test(vect<3>, vec_tup);
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

