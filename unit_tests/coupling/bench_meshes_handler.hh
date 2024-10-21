#ifndef BENCH_MESHES_HANDLER_HH_
#define BENCH_MESHES_HANDLER_HH_

#include <string>
#include <vector>
#include "input/reader_to_storage.hh"
#include "input/accessors.hh"
#include "input/input_type.hh"
#include "system/file_path.hh"


/**
 * Defines and provides path to mesh files used in assembly benchmark tests.
 */
class BenchMeshesHandler {
public:
    /// Declare Input::Type structure
    static Input::Type::Record & get_input_type() {
        std::string equation_name = "TestEquation";
        return Input::Type::Record("MeshesData", "List of benchmark meshes.")
            .declare_key("dg_asm",
                         Input::Type::Array( Input::Type::String(), 1, 12 ),
                         Input::Type::Default::obligatory(),
                         "Meshes of DG_asm test.")
            .declare_key("elasticity_asm",
                         Input::Type::Array( Input::Type::String(), 1, 12 ),
                         Input::Type::Default::obligatory(),
                         "Meshes of elasticity_asm test.")
            .declare_key("meshes_sizes",
                         Input::Type::Array(Input::Type::String(), 1, 3),
                         Input::Type::Default::obligatory(),
                         "Sizes of meshes.")
            .close();
    }

    /// Constructor
    BenchMeshesHandler() {
        std::string input_file = "../../../unit_tests/coupling/bench_meshes_handler.yaml";
        FilePath fp(input_file, FilePath::FileType::input_file);
        Input::ReaderToStorage reader( fp, BenchMeshesHandler::get_input_type() );
        in_rec_ = reader.get_root_interface<Input::Record>();
    }

    /// Return meshes of benchmark test given by asm_name
    std::vector<std::string> get_mesh_names(std::string asm_name) {
        Input::Array i_arr = in_rec_.val<Input::Array>(asm_name);
        std::vector<std::string> mesh_names;
    	for (Input::Iterator<std::string> it = i_arr.begin<std::string>();
                        it != i_arr.end();
                        ++it) {
    	    mesh_names.push_back(*it);
    	}
    	return mesh_names;
    }

    /// Return sizes of meshes
    std::vector<std::string> get_mesh_sizes() {
        Input::Array i_arr = in_rec_.val<Input::Array>("meshes_sizes");
        std::vector<std::string> mesh_sizes;
    	for (Input::Iterator<std::string> it = i_arr.begin<std::string>();
                        it != i_arr.end();
                        ++it) {
    	    mesh_sizes.push_back(*it);
    	}
    	return mesh_sizes;
    }

private:
    Input::Record in_rec_;
};


#endif /* BENCH_MESHES_HANDLER_HH_ */
