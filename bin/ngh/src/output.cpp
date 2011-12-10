#include "output.h"

#include "system.h"
#include "config.h"

void Output(char* output_fname, TMesh* mesh) {
    FILE* out = fopen(output_fname, "wt");
    if (out == NULL) {
        mythrow((char*) "Couldn't open output file.", __LINE__, __FUNC__);
    }

    fprintf(out, "$NeighbourFormat\n");
    fprintf(out, "1.0 0 %d\n", sizeof ( double));
    fprintf(out, "$EndNeighbourFormat\n");
    fprintf(out, "$Neighbours\n");
    fprintf(out, "%d\n", mesh->GetNNeighs());

    TNeighbour* ngh;

    FOR_NEIGHBOURS(ngh) {
        switch (ngh->GetType()) {
            case nt_bb:
                WriteNeighboursBB(ngh, out);
                break;
            case nt_vb:
                WriteNeighboursVB(ngh, out);
                break;
            case nt_vv:
                WriteNeighboursVV(ngh, out);
                break;
            default:
                mythrow((char*) "Runtime error - deny point.", __LINE__, __FUNC__);
        }
    }
    fprintf(out, "$EndNeighbours\n");
    fclose(out);
}
//=============================================================================
//
//=============================================================================

void WriteNeighboursBB(TNeighbour* ngh, FILE* out) {
    fprintf(out, "%d   %d ", ngh->GetId(), (int) ngh->GetType());
    fprintf(out, "%d ", ngh->GetNElements());

    int li;

    FOR_NEIGH_ELEMENTS(ngh, li) {
        fprintf(out, "%d %d ", ngh->getElement(li)->getLabel(), ngh->getSide(li)->GetLnum());
    }
    fprintf(out, "\n");
}
//=============================================================================
//
//=============================================================================

void WriteNeighboursVB(TNeighbour* ngh, FILE* out) {
    fprintf(out, "%d   %d ", ngh->GetId(), (int) ngh->GetType());
    //fprintf( out, "%d %d %d %lg\n", ngh->eid[ 0 ], ngh->eid[ 1 ], ngh->sid[ 1 ], ngh->GetCoef() );
    fprintf(out, "%d %d %d 1.0\n", ngh->getElement(0)->getLabel(), ngh->getElement(1)->getLabel(), ngh->getSide(1)->GetLnum());
    //fprintf(out, "%d %d %d %lg\n", ngh->getElement(0)->getLabel(), ngh->getElement(1)->getLabel(), ngh->getSide(1)->GetId(), ngh->GetCoef());
}

void WriteNeighboursVV(TNeighbour* ngh, FILE* out) {
    fprintf(out, "%d   %d ", ngh->GetId(), (int) ngh->GetType());
    fprintf(out, "%d %d %lg 1.00\n", ngh->getElement(0)->getLabel(), ngh->getElement(1)->getLabel(), ngh->GetCoef());
}

