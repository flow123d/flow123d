/*
 * element_data_cache_test.cpp
 *
 *  Created on: Jan 25, 2019
 *      Author: David Flanderka
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include "io/element_data_cache_base.hh"
#include "io/element_data_cache.hh"
#include "system/tokenizer.hh"
#include "la/distribution.hh"


/* Tests of ElementDataCacheBase functionality */
TEST(ElementDataCache, base_data_cache)
{
    // Default constructor and setters.
    {
        ElementDataCache<double> cache;
        cache.set_field_units("m*s^-2");
        cache.set_n_values(10);
        cache.set_dof_handler_hash(30);
        EXPECT_EQ(cache.field_units(), "m*s^-2");
        EXPECT_EQ(cache.n_values(), 10);
        EXPECT_EQ(cache.dof_handler_hash(), 30);
    }

    // 'Input' constructor.
    {
        ElementDataCache<double> cache("field_cache", 1.0, 1, 6);
        EXPECT_EQ(cache.field_input_name(), "field_cache");
        EXPECT_EQ(cache.get_time(), 1.0);
        EXPECT_TRUE(cache.is_actual(1.0, "field_cache"));
        EXPECT_FALSE(cache.is_actual(1.5, "field_cache"));
    }

    // 'Output' constructor.
    {
        ElementDataCache<double> cache("field_cache", 1, 3);
        EXPECT_EQ(cache.field_input_name(), "field_cache");
        EXPECT_EQ(cache.vtk_type(), ElementDataCacheBase::VTKValueType::VTK_FLOAT64);
        EXPECT_EQ(cache.n_values(), 3);
        EXPECT_EQ(cache.n_comp(), 1);
    }
}


TEST(ElementDataCache, read_data)
{
    std::stringstream ss; ss << "0 1 2 3 4 5 6 7 8 \n";
    Tokenizer tok(ss);
    ElementDataCache<double> data_cache("in_cache", 0.0, 1, 9);
    tok.next_line();
    for (unsigned int i=0; i<9; ++i) data_cache.read_ascii_data(tok, 1, i);

    auto &data_vec = *( data_cache.get_component_data(0).get() );
    for (unsigned int i=0; i<data_vec.size(); ++i) EXPECT_DOUBLE_EQ(data_vec[i], (double)i);

    EXPECT_EQ(data_cache.check_values(0.0, 0.0, 10.0), CheckResult::ok);

    data_cache.scale_data(0.1);
    for (unsigned int i=0; i<data_vec.size(); ++i) EXPECT_DOUBLE_EQ(data_vec[i], 0.1*i);
}


TEST(ElementDataCache, print_data)
{
	std::string expect_output = "1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 ";
	ElementDataCache<double> data_cache("out_cache", 1, 10);
	for (unsigned int i=0; i<data_cache.n_values(); ++i) {
		data_cache[i] = 1.0 + i*0.5;
    }

    double min, max;
    data_cache.get_min_max_range(min, max);
    EXPECT_DOUBLE_EQ(min, 1.0);
    EXPECT_DOUBLE_EQ(max, 5.5);

	{
        // print_ascii method
    	std::stringstream ss;
    	for (unsigned int i=0; i<data_cache.n_values(); ++i) data_cache.print_ascii(ss, i);
    	EXPECT_EQ(ss.str(), expect_output);
	}
	{
        // print_ascii_all method
    	std::stringstream ss;
    	data_cache.print_ascii_all(ss);
    	EXPECT_EQ(ss.str(), expect_output);
	}
}


TEST(ElementDataCache, value_operations)
{
	ElementDataCache<double> data_cache("data_cache", 3, 3);
	auto &data_vec = *( data_cache.get_component_data(0).get() );
	double val[] = { 0.0, 1.0, 2.0 };

	for (unsigned int i=0; i<data_cache.n_values(); ++i) data_cache.store_value(i, val);
	for (unsigned int i=0; i<data_vec.size(); ++i) EXPECT_DOUBLE_EQ( data_vec[i], (double)(i%3) );

	data_cache.add(2, val);
	data_cache.zero(0);
	for (unsigned int i=0; i<data_vec.size(); ++i) EXPECT_DOUBLE_EQ( data_vec[i], (double)(i%3)*(unsigned int)(i/3) );

	data_cache.add(0, val);
	data_cache.normalize(2, 2);
	for (unsigned int i=0; i<data_vec.size(); ++i) EXPECT_DOUBLE_EQ( data_vec[i], (double)(i%3) );
}


TEST(ElementDataCache, cache_size_operations)
{
	std::vector<unsigned int> offset_vec{ 0, 2, 5, 8, 12 };
	std::vector<unsigned int> opt_size_vec{ 0, 1, 0, 1, 2, 0, 1, 3, 0, 1, 2, 3 };
	std::vector<unsigned int> fix_size_vec{ 0, 1, 0, 0, 0, 1, 2, 0, 0, 1, 3, 0, 0, 1, 2, 3 };
	std::vector<double> avrg_vec{ 3.75, 4.75, 7, 9 };
	ElementDataCache<unsigned int> idx_cache("int_cache", 1, 12);
	for (unsigned int i=0; i<idx_cache.n_values(); ++i) idx_cache[i] = opt_size_vec[i];

	std::shared_ptr< ElementDataCache<unsigned int> > fixed_size_cache =
	        std::dynamic_pointer_cast< ElementDataCache<unsigned int> >(idx_cache.element_node_cache_fixed_size(offset_vec));
	for (unsigned int i=0; i<fixed_size_cache->n_values(); ++i) EXPECT_EQ( (*fixed_size_cache)[i], fix_size_vec[i] );

	std::shared_ptr< ElementDataCache<unsigned int> > optim_size_cache =
	        std::dynamic_pointer_cast< ElementDataCache<unsigned int> >(fixed_size_cache->element_node_cache_optimize_size(offset_vec));
	for (unsigned int i=0; i<optim_size_cache->n_values(); ++i) EXPECT_EQ( (*optim_size_cache)[i], opt_size_vec[i] );

	ElementDataCache<double> data_cache("double_cache", 1, 12);
	for (unsigned int i=0; i<data_cache.n_values(); ++i) data_cache[i] = (double)i;
	std::shared_ptr< ElementDataCache<double> > avrg_cache =
	        std::dynamic_pointer_cast< ElementDataCache<double> >(data_cache.compute_node_data(opt_size_vec, 4));
	for (unsigned int i=0; i<avrg_cache->n_values(); ++i) EXPECT_EQ( (*avrg_cache)[i], avrg_vec[i] );
}


TEST(ElementDataCache, gather)
{
	std::vector<double> vals_vec{ 0.5, 1.5, 2.0, 1.0, 0.8, 1.6, 2.4, 0.0, 2.5, 0.2, 3.7, 3.2 };
    Distribution * distr = new Distribution(4, MPI_COMM_WORLD);
    int rank = distr->myp();
    int n_proc = distr->np();

    LongIdx local_to_global[4];
    local_to_global[0]=2*rank; local_to_global[1]=2*rank+1; local_to_global[2]=rank+2*n_proc; local_to_global[3]=rank+3*n_proc;

    ElementDataCache<double> local_cache("field_cache", 1, 4);
    for (unsigned int i=0; i<local_cache.n_values(); ++i) local_cache[i] = vals_vec[ local_to_global[i] ];
    auto gathered_cache = local_cache.gather(distr, local_to_global);
    if (rank == 0) {
    	std::shared_ptr< ElementDataCache<double> > serial_cache = std::dynamic_pointer_cast< ElementDataCache<double> >(gathered_cache);
    	for (unsigned int i=0; i<serial_cache->n_values(); ++i) EXPECT_EQ( (*serial_cache)[i], vals_vec[i] );
    }
}
