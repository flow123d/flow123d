/*
 * reader_instances_test.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: jb
 */

#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <sstream>
#include <string>
#include <vector>
#include <mesh_constructor.hh>

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"
#include "io/reader_instances.hh"



Input::Record get_input_record(const std::string &input_str, Input::FileFormat format = Input::FileFormat::format_JSON) {

	istringstream is(input_str);
    Input::ReaderToStorage reader;
    IT::Record &in_rec = const_cast<IT::Record &>(Mesh::get_input_type());
    in_rec.finish();
    reader.read_stream(is, in_rec, format);

    return reader.get_root_interface<Input::Record>();
}


TEST(ReaderInstances, get_element_data) {
	Profiler::initialize();
	unsigned int i, j;

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Mesh * mesh = mesh_constructor();
    Input::Record i_rec = get_input_record("{mesh_file=\"fields/simplest_cube_data.msh\"}");
    auto reader = ReaderInstances::instance()->get_reader( i_rec );
    reader->read_physical_names(mesh);
    reader->read_raw_mesh(mesh);

    std::vector<int> el_ids;
    for (i=1; i<14; ++i) el_ids.push_back(i);


    // read data by components for MultiField
    bool actual_data = false;
    for (i=0; i<3; ++i) {
        typename ElementDataCache<int>::ComponentDataPtr multifield_data =
        		ReaderInstances::instance()->get_reader(i_rec)->get_element_data<int>("vector_fixed", 0.0, 13, 1,
        		actual_data, el_ids, i);
    	std::vector<int> &vec = *( multifield_data.get() );
    	EXPECT_EQ(13, vec.size());
    	for (j=0; j<mesh->element.size(); j++) EXPECT_EQ( i+1, vec[j] ); // bulk elements
    	for ( ; j<vec.size(); j++) EXPECT_EQ( i+4, vec[j] ); // boundary elements
    }


    // read data to one vector for Field
    actual_data=false;
    {
    	typename ElementDataCache<int>::ComponentDataPtr field_data =
        		ReaderInstances::instance()->get_reader(i_rec)->get_element_data<int>("vector_fixed", 1.0, 13, 3,
                actual_data, el_ids, 0);
    	std::vector<int> &vec = *( field_data.get() );
    	EXPECT_EQ(39, vec.size());
    	for (j=0; j<3*mesh->element.size(); j++) EXPECT_EQ( 2+(j%3), vec[j] ); // bulk elements
    	for ( ; j<vec.size(); j++) EXPECT_EQ( 5+(j%3), vec[j] ); // boundary elements
    }

    delete mesh;
}


TEST(ReaderInstances, find_header) {
    Profiler::initialize();

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Input::Record i_rec = get_input_record("{mesh_file=\"fields/simplest_cube_data.msh\"}");

    unsigned int n_elements=13;
    unsigned int n_comp=3;
    std::shared_ptr< std::vector<double> > data;
    std::vector<int> element_id_map(n_elements);
    for(unsigned int i=0;i<n_elements;i++)
    	element_id_map[i]=i+1;

    bool actual_data = false;
    data = ReaderInstances::instance()->get_reader(i_rec)->get_element_data<double>("vector_fixed", 0.0, n_elements, n_comp,
            actual_data, element_id_map, 0);
    EXPECT_EQ(1.0, (*data)[0]);
    EXPECT_EQ(2.0, (*data)[1]);
    EXPECT_EQ(3.0, (*data)[2]);
    EXPECT_EQ(3.0, (*data)[3*8+2]);
    EXPECT_EQ(4.0, (*data)[3*9]);

    data = ReaderInstances::instance()->get_reader(i_rec)->get_element_data<double>("vector_fixed", 0.1, n_elements, n_comp,
        actual_data, element_id_map, 0);
    EXPECT_EQ(1.0, (*data)[0]);
    EXPECT_EQ(2.0, (*data)[1]);
    EXPECT_EQ(3.0, (*data)[2]);
    EXPECT_EQ(3.0, (*data)[3*8+2]);
    EXPECT_EQ(4.0, (*data)[3*9]);

    data = ReaderInstances::instance()->get_reader(i_rec)->get_element_data<double>("vector_fixed", 0.9, n_elements, n_comp,
            actual_data, element_id_map, 0);
    EXPECT_EQ(1.0, (*data)[0]);
    EXPECT_EQ(2.0, (*data)[1]);
    EXPECT_EQ(3.0, (*data)[2]);
    EXPECT_EQ(3.0, (*data)[3*8+2]);
    EXPECT_EQ(4.0, (*data)[3*9]);

    data = ReaderInstances::instance()->get_reader(i_rec)->get_element_data<double>("vector_fixed", 1.0, n_elements, n_comp,
            actual_data, element_id_map, 0);
    EXPECT_EQ(2.0, (*data)[0]);
    EXPECT_EQ(3.0, (*data)[1]);
    EXPECT_EQ(4.0, (*data)[2]);
    EXPECT_EQ(4.0, (*data)[3*8+2]);
    EXPECT_EQ(5.0, (*data)[3*9]);

    data = ReaderInstances::instance()->get_reader(i_rec)->get_element_data<double>("vector_fixed", 1.1, n_elements, n_comp,
            actual_data, element_id_map, 0);
    EXPECT_EQ(2.0, (*data)[0]);
    EXPECT_EQ(3.0, (*data)[1]);
    EXPECT_EQ(4.0, (*data)[2]);
    EXPECT_EQ(4.0, (*data)[3*8+2]);
    EXPECT_EQ(5.0, (*data)[3*9]);

    data = ReaderInstances::instance()->get_reader(i_rec)->get_element_data<double>("vector_fixed", 2.1, n_elements, n_comp,
            actual_data, element_id_map, 0);
    EXPECT_EQ(2.0, (*data)[0]);
    EXPECT_EQ(3.0, (*data)[1]);
    EXPECT_EQ(4.0, (*data)[2]);
    EXPECT_EQ(4.0, (*data)[3*8+2]);
    EXPECT_EQ(5.0, (*data)[3*9]);

    data = ReaderInstances::instance()->get_reader(i_rec)->get_element_data<double>("vector_fixed", 200, n_elements, n_comp,
            actual_data, element_id_map, 0);
    EXPECT_EQ(2.0, (*data)[0]);
    EXPECT_EQ(3.0, (*data)[1]);
    EXPECT_EQ(4.0, (*data)[2]);
    EXPECT_EQ(4.0, (*data)[3*8+2]);
    EXPECT_EQ(5.0, (*data)[3*9]);

}


TEST(ReaderInstances, get_reader) {
	Profiler::initialize();

	// has to introduce some flag for passing absolute path to 'test_units' in source tree
	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

	{
		FilePath mesh1("mesh/test_input.msh", FilePath::input_file);
		FilePath mesh2("mesh/simplest_cube.msh", FilePath::input_file);
		Input::Record i_rec1 = get_input_record("{mesh_file=\"mesh/test_input.msh\"}");
		Input::Record i_rec2 = get_input_record("{mesh_file=\"mesh/simplest_cube.msh\"}");
		ReaderInstances::instance()->get_reader(i_rec1);
		ReaderInstances::instance()->get_reader(i_rec2);
	}

	{
	    Mesh * mesh = mesh_constructor();
		Input::Record i_rec = get_input_record("{mesh_file=\"mesh/test_input.msh\"}");
		auto mesh_reader = ReaderInstances::instance()->get_reader(i_rec);
		mesh_reader->read_physical_names(mesh);
		mesh_reader->read_raw_mesh(mesh);

		EXPECT_EQ(118, mesh->n_nodes());

		delete mesh;
	}

	{
	    Mesh * mesh = mesh_constructor();
		Input::Record i_rec = get_input_record("{mesh_file=\"mesh/simplest_cube.msh\"}");
		auto mesh_reader = ReaderInstances::instance()->get_reader(i_rec);
		mesh_reader->read_physical_names(mesh);
		mesh_reader->read_raw_mesh(mesh);

		EXPECT_EQ(8, mesh->n_nodes());

		delete mesh;
	}
}


TEST(ReaderInstances, repeat_call) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Input::Record i_rec = get_input_record("{mesh_file=\"mesh/test_108_elem.msh\"}");
    for (unsigned int i=0; i<2; ++i) {
        auto mesh_reader = ReaderInstances::instance()->get_reader(i_rec);
        Mesh * mesh = mesh_constructor();
        mesh_reader->read_physical_names(mesh);
        mesh_reader->read_raw_mesh(mesh);
        delete mesh;
    }
}
