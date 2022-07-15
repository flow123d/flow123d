/*
 * gmsh_reader_test.cpp
 *
 *  Created on: Oct 2, 2012
 *      Author: jb
 */

#define FEAL_OVERRIDE_ASSERTS

#include <flow_gtest.hh>
#include <sstream>
#include <string>
#include <mesh_constructor.hh>
#include <gmsh.h>

#include "system/sys_profiler.hh"
#include "system/file_path.hh"

#include "mesh/mesh.h"
#include "io/msh_gmshreader.h"

class Gmsh4MeshReader
{
public:
    Gmsh4MeshReader(const FilePath &mesh_file) {
        gmsh::initialize();
        gmsh::open( mesh_file.filename() );
    }

    ~Gmsh4MeshReader() {
        gmsh::clear();
        gmsh::finalize();
    }

    void read_physical_names(Mesh * mesh) {
        ASSERT_PTR(mesh).error("Argument mesh is NULL.\n");

        gmsh::vectorpair dim_tags;
        gmsh::model::getPhysicalGroups(dim_tags);

        for (auto dim_tag : dim_tags) {
            std::string physical_name;
            gmsh::model::getPhysicalName(dim_tag.first, dim_tag.second, physical_name);
            mesh->add_physical_name( dim_tag.first, dim_tag.second, physical_name );
        }
    }

    void read_nodes(Mesh * mesh) {
        ASSERT_PTR(mesh).error("Argument mesh is NULL.\n");

        MessageOut() << "- Reading nodes...";

        std::vector<std::size_t> node_tags;
        std::vector<double> coord;
        std::vector<double> parametric_coord;
        gmsh::model::mesh::getNodes(node_tags, coord, parametric_coord);

        unsigned int n_nodes = node_tags.size();
        mesh->init_node_vector( n_nodes );
        //if (n_nodes == 0) THROW( ExcZeroNodes() << EI_Position(tok_.position_msg()) ); // TODO fix exception

        for (uint i=0; i<n_nodes; i++) {
            unsigned int node_id = (unsigned int)(node_tags[i]);
            unsigned int offset = 3*i;
            arma::vec3 coords;
            coords(0) = coord[offset];
            coords(1) = coord[offset + 1];
            coords(2) = coord[offset + 2];
            mesh->add_node(node_id, coords);
        }

        MessageOut().fmt("... {} nodes read. \n", n_nodes);
    }

    void read_elements(Mesh * mesh) {
        ASSERT_PTR(mesh).error("Argument mesh is NULL.\n");

        MessageOut() << "- Reading elements...";

        gmsh::vectorpair dim_tags;
        gmsh::model::getEntities(dim_tags);

        std::vector<int> element_types;
        std::vector<std::vector<std::size_t> > element_tags;
        std::vector<std::vector<std::size_t> > node_tags;
        gmsh::model::mesh::getElements(element_types, element_tags, node_tags);
        uint n_elements = 0;
        for (auto t : element_tags) n_elements += t.size();
        mesh->init_element_vector(n_elements);

        for (auto dim_tag : dim_tags) {
            int entity_dim = dim_tag.first;
            int entity_tag = dim_tag.second;
            std::vector<int> physical_tags;
            gmsh::model::getPhysicalGroupsForEntity(entity_dim, entity_tag, physical_tags);

            if (physical_tags.size() == 0) continue; // no elements in entity

            int region_id = physical_tags[0];
            gmsh::model::mesh::getElements(element_types, element_tags, node_tags, entity_dim, entity_tag);

            for (uint i=0; i<element_types.size(); ++i) {
                unsigned int dim;
                switch (element_types[i]) {
                    case 1:
                        dim = 1;
                        break;
                    case 2:
                        dim = 2;
                        break;
                    case 4:
                        dim = 3;
                        break;
                    case 15:
                        dim = 0;
                        break;
                    default:
                        dim = 0;
                        // TODO fix throw, set correct EI_GMSHFile
                        //THROW(ExcUnsupportedType() << EI_ElementId(id) << EI_ElementType(type) << EI_GMSHFile(tok_.f_name()) );
                        break;
                }
                unsigned int partition_id = 0; // not supported in Gmsh4
                std::vector<unsigned int> node_ids(dim+1);
                for (uint j=0; j<element_tags[i].size(); ++j) {
                    unsigned int id = (unsigned int)(element_tags[i][j]);
                    unsigned int offset = (dim+1) * j;
                    for (uint n=0; n<dim+1; ++n)
                        node_ids[n] = node_tags[i][offset+n];
                    mesh->add_element(id, dim, region_id, partition_id, node_ids);
                }
            }
        }

        MessageOut().fmt("... {} bulk elements, {} boundary elements. \n", mesh->n_elements(), mesh->bc_mesh()->n_elements());
    }
};


TEST(GMSHReader, read_mesh_from_file) {
    Profiler::instance();

    // has to introduce some flag for passing absolute path to 'test_units' in source tree
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    std::string mesh_in_string = "{mesh_file=\"mesh/square_2frac_msh4.msh\"}";
    FilePath mesh_file("mesh/square_2frac_msh4.msh", FilePath::FileType::input_file);
    Mesh * mesh = mesh_constructor(mesh_in_string);
    Gmsh4MeshReader reader(mesh_file);
    reader.read_physical_names(mesh);
//    reader->read_raw_mesh(mesh);
    reader.read_nodes(mesh);
    reader.read_elements(mesh);
//    mesh->setup_topology();
//    mesh->check_and_finish();

//    EXPECT_EQ(9, mesh->region_db().size());
//    EXPECT_EQ(3, mesh->region_db().bulk_size());
    EXPECT_EQ(151, mesh->n_nodes());
    EXPECT_EQ(270, mesh->n_elements()); // 270 bulk + 42 boundary elements = 312 elements total

    delete mesh;
}
