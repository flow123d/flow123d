//---------------------------------------------------------------------------
#include <iostream>

#include "config.h"
#include "problem.h"
#include "output.h"
#include "TInstanceProfiler.h"
#include "TMeshReader.h"
#include "TNeighboursVV.h"
//---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    std::cout << "This is NGH generator, version " << _VERSION_ << "\n";
    std::cout << "Built on " << __DATE__ << " at " << __TIME__ << ".\n\n";

    if (argc < 2) {
        std::cout << "Usage: ngh.exe ini_file\n";
        return 0;
    }

    /** Zrizeni instance TMesh */
    TMesh* mesh = new TMesh();
    /** Zrizeni instance TMeshReader */
    TMeshReader* meshReader = new TMeshReader();

    try {
        TProblem* problem = new TProblem(argv[ 1 ]);

        meshReader->read(problem->getMesh_fname(), mesh);

        mesh->CreateSides();

        mesh->NodeToElement();
        mesh->SideToNode();

        mesh->CreateEdges();

        std::cout << "Searching for neighbours - Start\n";
        if (problem->getBB()) mesh->CreateNeighboursBB();
        if (problem->getVB()) mesh->CreateNeighboursVB();
        if (problem->getVV()) {
            // mesh->CreateGeometry();
            // mesh->TestElements();
            TNeighboursVV* neigVV = new TNeighboursVV();
            neigVV->create(problem, mesh);
        }
        std::cout << "Searching for neighbours - End\n";

        Output(problem->getOutput_fname(), mesh);
        std::cout << "\nCalculation finished succesfuly.\nI'm going home.\n";
        if (problem->getPause_after_run()) {
            std::cin.get();
        }

        delete problem;
    } catch (char *ch) {
        std::cout << ch << "\n";
        return 1;
    }

    TInstanceProfiler::printResult();

    delete mesh;
    delete meshReader;

    return 0;
}
