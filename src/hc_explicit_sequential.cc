#include "hc_explicit_sequential.hh"
#include "flow/darcy_flow_mh.hh"
#include "flow/darcy_flow_mh_output.hh"
#include "transport_operator_splitting.hh"
#include "transport.h"
#include "equation.hh"
#include "time_marks.hh"
#include "mesh/mesh.h"
#include "mesh/msh_gmshreader.h"
#include "io/output.h"

/**
 * FUNCTION "MAIN" FOR COMPUTING MIXED-HYBRID PROBLEM FOR UNSTEADY SATURATED FLOW
 */
HC_ExplicitSequential::HC_ExplicitSequential(ProblemType problem_type)
{
    type_=problem_type;

    // Initialize Time Marks
    main_time_marks = new TimeMarks();

    // Material Database
    const string& material_file_name = IONameHandler::get_instance()->get_input_file_name(OptGetStr( "Input", "Material", "\\" ));
    material_database = new MaterialDatabase(material_file_name);

    // Read mesh
    mesh = new Mesh();
    const string& mesh_file_name = IONameHandler::get_instance()->get_input_file_name(OptGetStr("Input", "Mesh", NULL));
    MeshReader* meshReader = new GmshMeshReader();
    meshReader->read(mesh_file_name, mesh);
    mesh->setup_topology();
    mesh->setup_materials(*material_database);
    Profiler::instance()->set_task_size(mesh->n_elements());

    // setup water flow object
    switch (problem_type) {
        case STEADY_SATURATED:
            water=new DarcyFlowMH_Steady(main_time_marks, mesh, material_database);
            break;
        case UNSTEADY_SATURATED:
            water=new DarcyFlowMH_Unsteady(main_time_marks, mesh, material_database);
            break;
        case UNSTEADY_SATURATED_LMH:
            water=new DarcyFlowLMH_Unsteady(main_time_marks, mesh, material_database);
            break;
        default:
            xprintf(Err,"Wrong problem type: %d",problem_type);
            break;
    }
    // object for water postprocessing and output
    water_output = new DarcyFlowMHOutput(water);

    // optionally setup transport objects
    if ( OptGetBool("Transport", "Transport_on", "no") ) {
        transport_reaction = new TransportOperatorSplitting(material_database, mesh);
    } else {
        transport_reaction = new EquationNothing();
    }



}

void HC_ExplicitSequential::run_simulation()
{
    OutputTime *output_time;

    int i, rank;
    // setup output
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if(rank == 0) {
        string output_file = IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Output", "Output_file", "\\"));
        output_time = new OutputTime(mesh, output_file);
    }

    // ensure we have planned times for both processes
    water->choose_next_time();
    transport_reaction->choose_next_time();
    while (! (water->is_end() && transport_reaction->is_end() ) ) {

        //
        water->compute_one_step();
        water_output->postprocess();
        water_output->output();

    }
}


//-----------------------------------------------------------------------------
// vim: set cindent:
//-----------------------------------------------------------------------------
