/*
 * gmsh_reader_test.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: jb
 */


#include <flow_gtest.hh>
#include <sstream>
#include <string>
#include <vector>

#include "system/sys_profiler.hh"

#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "mesh/reader_instances.hh"



TEST(GMSHReader, read_mesh_from_stream) {
    Profiler::initialize();
    
    string fname = string(UNIT_TESTS_SRC_DIR) + "/mesh/simplest_cube.msh";
    ifstream ifs( fname.c_str() );
    if (! ifs) cout << "Can not open file!" << endl;
    stringstream ss;
    ss << ifs.rdbuf();

    Mesh mesh;
    GmshMeshReader reader(ss);

    reader.read_mesh(&mesh);

    EXPECT_EQ(9, mesh.n_elements());

}


TEST(GMSHReader, read_mesh_from_file) {
    Profiler::initialize();

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    FilePath mesh_file("mesh/test_input.msh", FilePath::input_file);


    Mesh mesh;
    GmshMeshReader reader(mesh_file);

    reader.read_mesh(&mesh);

    EXPECT_EQ(216, mesh.n_elements());

}


TEST(GMSHReader, get_element_data) {
	Profiler::initialize();
	unsigned int i, j;

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    FilePath mesh_file("fields/simplest_cube_data.msh", FilePath::input_file);

    Mesh mesh;
    ReaderInstances::instance()->get_reader(mesh_file)->read_mesh(&mesh);

    std::vector<int> el_ids;
    for (i=1; i<14; ++i) el_ids.push_back(i);


    // read data by components for MultiField
    GMSH_DataHeader search_header;
    search_header.actual=false;
    search_header.field_name="vector_fixed";
    search_header.n_components=1;
    search_header.n_entities=13;
    search_header.time=0.0;

    for (i=0; i<3; ++i) {
        typename ElementDataCache<int>::ComponentDataPtr multifield_data =
        		ReaderInstances::instance()->get_reader(mesh_file)->get_element_data<int>(search_header, el_ids, i);
    	std::vector<int> &vec = *( multifield_data.get() );
    	EXPECT_EQ(13, vec.size());
    	for (j=0; j<mesh.element.size(); j++) EXPECT_EQ( i+1, vec[j] ); // bulk elements
    	for ( ; j<vec.size(); j++) EXPECT_EQ( i+4, vec[j] ); // boundary elements
    }


    // read data to one vector for Field
    search_header.actual=false;
    search_header.n_components=3;
    search_header.time=1.0;

    {
    	typename ElementDataCache<int>::ComponentDataPtr field_data =
        		ReaderInstances::instance()->get_reader(mesh_file)->get_element_data<int>(search_header, el_ids, 0);
    	std::vector<int> &vec = *( field_data.get() );
    	EXPECT_EQ(39, vec.size());
    	for (j=0; j<3*mesh.element.size(); j++) EXPECT_EQ( 2+(j%3), vec[j] ); // bulk elements
    	for ( ; j<vec.size(); j++) EXPECT_EQ( 5+(j%3), vec[j] ); // boundary elements
    }
}


TEST(GMSHReader, find_header) {
    Profiler::initialize();

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    FilePath mesh_file("fields/simplest_cube_data.msh", FilePath::input_file);

    GmshMeshReader reader(mesh_file);
    unsigned int n_elements=13;
    unsigned int n_comp=3;
    std::shared_ptr< std::vector<double> > data;
    std::vector<int> element_id_map(n_elements);
    for(unsigned int i=0;i<n_elements;i++)
    	element_id_map[i]=i+1;
    GMSH_DataHeader header;
    header.actual=false;
    header.field_name="vector_fixed";
    header.n_components=n_comp;
    header.n_entities=n_elements;

    header.time=0.0;
    data = reader.get_element_data<double>(header, element_id_map, 0);
    EXPECT_EQ(1.0, (*data)[0]);
    EXPECT_EQ(2.0, (*data)[1]);
    EXPECT_EQ(3.0, (*data)[2]);
    EXPECT_EQ(3.0, (*data)[3*8+2]);
    EXPECT_EQ(4.0, (*data)[3*9]);

    header.time=0.1;
    data = reader.get_element_data<double>(header, element_id_map, 0);
    EXPECT_EQ(1.0, (*data)[0]);
    EXPECT_EQ(2.0, (*data)[1]);
    EXPECT_EQ(3.0, (*data)[2]);
    EXPECT_EQ(3.0, (*data)[3*8+2]);
    EXPECT_EQ(4.0, (*data)[3*9]);

    header.time=0.9;
    data = reader.get_element_data<double>(header, element_id_map, 0);
    EXPECT_EQ(1.0, (*data)[0]);
    EXPECT_EQ(2.0, (*data)[1]);
    EXPECT_EQ(3.0, (*data)[2]);
    EXPECT_EQ(3.0, (*data)[3*8+2]);
    EXPECT_EQ(4.0, (*data)[3*9]);

    header.time=1.0;
    data = reader.get_element_data<double>(header, element_id_map, 0);
    EXPECT_EQ(2.0, (*data)[0]);
    EXPECT_EQ(3.0, (*data)[1]);
    EXPECT_EQ(4.0, (*data)[2]);
    EXPECT_EQ(4.0, (*data)[3*8+2]);
    EXPECT_EQ(5.0, (*data)[3*9]);

    header.time=1.1;
    data = reader.get_element_data<double>(header, element_id_map, 0);
    EXPECT_EQ(2.0, (*data)[0]);
    EXPECT_EQ(3.0, (*data)[1]);
    EXPECT_EQ(4.0, (*data)[2]);
    EXPECT_EQ(4.0, (*data)[3*8+2]);
    EXPECT_EQ(5.0, (*data)[3*9]);

    header.time=2.1;
    data = reader.get_element_data<double>(header, element_id_map, 0);
    EXPECT_EQ(2.0, (*data)[0]);
    EXPECT_EQ(3.0, (*data)[1]);
    EXPECT_EQ(4.0, (*data)[2]);
    EXPECT_EQ(4.0, (*data)[3*8+2]);
    EXPECT_EQ(5.0, (*data)[3*9]);

    header.time=200;
    data = reader.get_element_data<double>(header, element_id_map, 0);
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
		ReaderInstances::instance()->get_reader(mesh1);
		ReaderInstances::instance()->get_reader(mesh2);
	}

	{
	    Mesh mesh;
		FilePath mesh_file("mesh/test_input.msh", FilePath::input_file);
		std::shared_ptr<GmshMeshReader> mesh_reader = ReaderInstances::instance()->get_reader(mesh_file);

		mesh_reader->read_mesh(&mesh);
		EXPECT_EQ(118, mesh.n_nodes());
	}

	{
	    Mesh mesh;
		FilePath mesh_file("mesh/simplest_cube.msh", FilePath::input_file);
		std::shared_ptr<GmshMeshReader> mesh_reader = ReaderInstances::instance()->get_reader(mesh_file);

		mesh_reader->read_mesh(&mesh);
		EXPECT_EQ(8, mesh.n_nodes());
	}
}
