/*!
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file   postprocess.cc
 * @brief  Postprocessing
 *
 */

#include "constantdb.h"
#include "mesh/ini_constants_mesh.hh"

#include "system.hh"
#include "math_fce.h"
#include "darcy_flow_mh.hh"
#include "problem.h"
#include "mesh.h"
#include "boundaries.h"
#include "sources.h"
#include "transport.h"
#include "output.h"
#include "postprocess.h"
#include "materials.hh"
#include "transport_bcd.h"
#include "concentrations.h"

static void calc_external_balance(Mesh*);
static void make_side_flux(struct Problem*, Mesh*);
static void make_element_scalar(struct Problem*, Mesh*);
static void make_element_vector(struct Problem*, Mesh*);
static void make_element_vector_line(ElementFullIter);
static void make_element_vector_triangle(ElementFullIter);
static void make_element_vector_tetrahedron(ElementFullIter);
static void make_sides_scalar(struct Problem*, Mesh*);
static void make_node_scalar(Mesh*);
static void make_neighbour_flux(struct Problem*, Mesh*);
static void make_previous_scalar(struct Problem*, Mesh*);
static void calc_element_balance(Mesh*);
static void water_balance(Mesh*, MaterialDatabase*);
static double calc_water_balance(Mesh*, int);

//void make_node_vector(Mesh*);

//=============================================================================
// CONVERT SOLUTION, CALCULATE BALANCES, ETC...
//=============================================================================

void postprocess(struct Problem *problem) {
    ASSERT(!(problem == NULL), "NULL as argument of function postprocess()\n");
    xprintf(Msg, "Converting solution... ")/*orig verb 2*/;

    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    make_side_flux(problem, mesh);
    make_element_scalar(problem, mesh);
    make_element_vector(problem, mesh);
    make_sides_scalar(problem, mesh);
    make_node_scalar(mesh);
    //make_node_vector( mesh );
    make_neighbour_flux(problem, mesh);
    make_previous_scalar(problem, mesh);
    xprintf(Msg, "O.K.\n")/*orig verb 2*/;
    calc_element_balance(mesh);
    water_balance(mesh, problem->material_database);
    //  if (problem->transport_on == true)
    //         transport( problem );
    xprintf(Msg, "Postprocessing phase O.K.\n")/*orig verb 2*/;
}
//=============================================================================
// FILL TH "FLUX" FIELD FOR ALL SIDES IN THE MESH
//=============================================================================

void make_side_flux(struct Problem *problem, Mesh* mesh) {
    int li, soi;
    double *sol;
    ElementIter ele;
    struct Side *sde;

    soi = 0;
    sol = problem->water->solution_vector();
    FOR_ELEMENTS(ele)
    for (li = 0; li < ele->n_sides; li++) {
        sde = ele->side[ li ];
        sde->flux = sol[ soi++ ];
        //if( fabs( sde->flux ) < ZERO )
        //	sde->flux = 0.0;
    }
}
//=============================================================================
// FILL TH "SCALAR" FIELD FOR ALL ELEMENTS IN THE MESH
//=============================================================================

void make_element_scalar(struct Problem *problem, Mesh* mesh) {
    int soi;
    double *sol;
    ElementIter ele;

    soi = mesh->n_sides;
    sol = problem->water->solution_vector();
    FOR_ELEMENTS(ele)
    ele->scalar = sol[ soi++ ];
}

/****
 * compute Darcian velocity in centre of elements
 *
 */
void make_element_vector(struct Problem *problem, Mesh* mesh) {
    ElementIter ele;
    //FILE *out;

    // upravy Ji -- pomocny tisk
    //out = xfopen( "pomout2.txt", "wt" );
    //xfprintf( out, "Pomocny tisk bazovych funkci po vypoctu\n\n");

    FOR_ELEMENTS(ele) {
        switch (ele->type) {
            case LINE:
                make_element_vector_line(ele);
                break;
            case TRIANGLE:
                make_element_vector_triangle(ele);
                //xfprintf( out, "%d \n", ele->id);
                //xfprintf( out, "Plocha %12.8f\n", ele->metrics);
                //xfprintf( out, "Teziste\n");
                //xfprintf( out, "%12.8f %12.8f %12.8f \n", ele->centre[0], ele->centre[1], ele->centre[2]);
                //xfprintf( out, "Bazove funkce\n");
                //xfprintf( out, "%12.8f %12.8f %12.8f \n", ele->bas_alfa[0], ele->bas_beta[0], ele->bas_gama[0] );
                //xfprintf( out, "%12.8f %12.8f %12.8f \n", ele->bas_alfa[1], ele->bas_beta[1], ele->bas_gama[1] );
                //xfprintf( out, "%12.8f %12.8f %12.8f \n", ele->bas_alfa[2], ele->bas_beta[2], ele->bas_gama[2] );
                //xfprintf( out, " \n");
                //xfprintf( out, "Pretoky\n");
                //xfprintf( out, "%12.8f %12.8f %12.8f \n", ele->side[0]->flux, ele->side[1]->flux, ele->side[2]->flux);
                //xfprintf( out, "Vektor v tezisti\n");
                //xfprintf( out, "%12.8f %12.8f %12.8f \n",  ele->vector[0], ele->vector[1], ele->vector[2]);
                //xfprintf( out, " \n");
                //xfprintf( out, " \n\n");
                break;
            case TETRAHEDRON:
                make_element_vector_tetrahedron(ele);
                break;
        }
        ele->v_length = vector_length(ele->vector);
    }
    //xfclose( out );
}
//=============================================================================
//
//=============================================================================

void make_element_vector_line(ElementFullIter ele) {
    double darcy_vel = (ele->side[1]->flux - ele->side[0]->flux) / 2.0 / ele->material->size;

    // normalize element vector [node 0, node 1]
    ele->vector[ 0 ] = (ele->node[1]->getX() - ele->node[0]->getX());
    ele->vector[ 1 ] = (ele->node[1]->getY() - ele->node[0]->getY());
    ele->vector[ 2 ] = (ele->node[1]->getZ() - ele->node[0]->getZ());
    normalize_vector(ele->vector);

    ele->vector[ 0 ] *= darcy_vel;
    ele->vector[ 1 ] *= darcy_vel;
    ele->vector[ 2 ] *= darcy_vel;

}
//=============================================================================
//
//=============================================================================

void make_element_vector_triangle(ElementFullIter ele) {
    double bas[ 3 ][ 3 ];
    double ex[ 3 ];
    double ey[ 3 ];
    double ez[ 3 ];
    double tmp[ 3 ];
    double X[ 3 ];
    double mid[3];
    double u[3];
    int i, li;

    // begin -- upravy Ji. -- prepocet vektoru do teziste elementu
    // 23.2.2007
    ex[ 0 ] = ele->node[1]->getX() - ele->node[0]->getX();
    ex[ 1 ] = ele->node[1]->getY() - ele->node[0]->getY();
    ex[ 2 ] = ele->node[1]->getZ() - ele->node[0]->getZ();

    tmp[ 0 ] = ele->node[2]->getX() - ele->node[0]->getX();
    tmp[ 1 ] = ele->node[2]->getY() - ele->node[0]->getY();
    tmp[ 2 ] = ele->node[2]->getZ() - ele->node[0]->getZ();
    vector_product(ex, tmp, ez);
    normalize_vector(ex);
    normalize_vector(ez);
    vector_product(ez, ex, ey);
    normalize_vector(ey);

    u[0] = ele->centre[0] - ele->node[0]->getX();
    u[1] = ele->centre[1] - ele->node[0]->getY();
    u[2] = ele->centre[2] - ele->node[0]->getZ();
    mid[0] = scalar_product(ex, u);
    mid[1] = scalar_product(ey, u);
    mid[2] = scalar_product(ez, u);


    for (li = 0; li < 3; li++)
        ele->vector[ li ] = 0;
    for (i = 0; i < 3; i++) {
        bas[ i ][ 0 ] = ele->bas_gama[ i ] * (mid[0] - ele->bas_alfa[ i ]);
        bas[ i ][ 1 ] = ele->bas_gama[ i ] * (mid[1] - ele->bas_beta[ i ]);
        bas[ i ][ 2 ] = 0;
        for (li = 0; li < 2; li++)
            ele->vector[ li ] += ele->side[ i ]->flux * bas[ i ][ li ] / ele->material->size;
    }

    /*for (li = 0; li < 3; li++){
            A[ 0 ][ li ] = ex[ li ];
            A[ 1 ][ li ] = ey[ li ];
            A[ 2 ][ li ] = ez[ li ];
    }
    matrix_x_matrix(ele->vector, 1, 3, A, 3, 3, X);
    for (li = 0; li < 3; li++)
            ele->vector[ li ] = X[ li ];

     */
    X[0] = ele->vector[0] * ex[0] +
            ele->vector[1] * ey[0] +
            ele->vector[2] * ez[0];
    X[1] = ele->vector[0] * ex[1] +
            ele->vector[1] * ey[1] +
            ele->vector[2] * ez[1];
    X[2] = ele->vector[0] * ex[2] +
            ele->vector[1] * ey[2] +
            ele->vector[2] * ez[2];
    for (li = 0; li < 3; li++)
        ele->vector[ li ] = X[ li ];
    // end -- upravy Ji. -- prepocet vektoru do teziste elementu

}
//=============================================================================
//
//=============================================================================

void make_element_vector_tetrahedron(ElementFullIter ele) {
    double bas[ 4 ][ 3 ];
    int i, li;

    for (li = 0; li < 3; li++)
        ele->vector[ li ] = 0;
    for (i = 0; i < 4; i++) {
        bas[ i ][ 0 ] = ele->bas_delta[ i ] * (ele->centre[ 0 ]
                - ele->bas_alfa[ i ]);
        bas[ i ][ 1 ] = ele->bas_delta[ i ] * (ele->centre[ 1 ]
                - ele->bas_beta[ i ]);
        bas[ i ][ 2 ] = ele->bas_delta[ i ] * (ele->centre[ 2 ]
                - ele->bas_gama[ i ]);
        for (li = 0; li < 3; li++)
            ele->vector[ li ] += ele->side[ i ]->flux * bas[ i ][ li ];
    }
}
//=============================================================================
// FILL TH "SCALAR" FIELD FOR ALL INTERNAL SIDES IN THE MESH
//=============================================================================

void make_sides_scalar(struct Problem *problem, Mesh* mesh) {
    struct Edge *edg;
    double *sol;
    int soi, si;
    struct Side *sde;

    soi = mesh->n_sides + mesh->n_elements();
    sol = problem->water->solution_vector();

    FOR_EDGES(edg) {
        for (si = 0; si < edg->n_sides; si++) {
            sde = edg->side[ si ];
            sde->scalar = sol[ soi ];
        }
        soi++;
    }
}
//=============================================================================
// pressure interpolation
//
// some sort of weighted average over elements and sides/edges of an node
//
// TODO:
// questions:
// - explain details of averaging and motivation for this type
// - edge scalars are stronger since thay are computed twice by each side
// - why division by (nod.aux-1)
// - ?implement without TED, TSD global arrays with single loop over elements, sides and one loop over nodes for normalization
//
//=============================================================================

void make_node_scalar(Mesh* mesh) {
    int li;
    NodeIter nod;
    struct Side *sde;
    double dist;
    int max_side_id = 0;

    double **TED;
    double **TSD;


    TED = (double **) xmalloc((mesh->element.size() + 1) * sizeof (double *));

    FOR_SIDES(sde)
    if (max_side_id <= sde->id)
        max_side_id = sde->id;

    TSD = (double **) xmalloc((max_side_id + 1) * sizeof (double *));

    FOR_ELEMENTS(ele)
    TED[ele.index()] = (double*) xmalloc(ele->n_nodes * sizeof (double));
    FOR_SIDES(sde)
    TSD[sde->id] = (double*) xmalloc(sde->n_nodes * sizeof (double));

    FOR_NODES( nod ) {
        nod->scalar = 0.0;
        nod->faux = 0.0;
        nod->aux = 0;
    }
    FOR_ELEMENTS(ele)
    for (li = 0; li < ele->n_nodes; li++) {
        nod = ele->node[li];

        dist = sqrt(
                ((nod->getX() - ele->centre[ 0 ])*(nod->getX() - ele->centre[ 0 ])) +
                ((nod->getY() - ele->centre[ 1 ])*(nod->getY() - ele->centre[ 1 ])) +
                ((nod->getZ() - ele->centre[ 2 ])*(nod->getZ() - ele->centre[ 2 ]))
                );

        TED[ele.index()][li] = dist;
        nod->faux += dist; //       nod->faux += 1 / dist;
        nod->aux++;
    }

    FOR_SIDES(sde) {
        for (li = 0; li < sde->n_nodes; li++) {
            nod = sde->node[li];

        dist = sqrt(
                ((nod->getX() - sde->centre[ 0 ])*(nod->getX() - sde->centre[ 0 ])) +
                ((nod->getY() - sde->centre[ 1 ])*(nod->getY() - sde->centre[ 1 ])) +
                ((nod->getZ() - sde->centre[ 2 ])*(nod->getZ() - sde->centre[ 2 ]))
                );

            TSD[sde->id][li] = dist;
            nod->faux += dist; //      nod->faux += 1 / dist;
            nod->aux++;
        }
    }
    FOR_ELEMENTS(ele)
    for (li = 0; li < ele->n_nodes; li++) {
        nod = ele->node[li];
        nod->scalar += ele->scalar * (1 - TED[ele.index()][li] / nod->faux)
                / (nod->aux - 1); // 1 / (dist * nod->faux);
    }
    FOR_SIDES(sde)
    for (li = 0; li < sde->n_nodes; li++) {
        nod = sde->node[li];
        nod->scalar += sde->scalar * (1 - TSD[sde->id][li] / nod->faux)
                / (nod->aux - 1); // 1 / (dist * nod->faux);
    }
    xfree(TED);
    xfree(TSD);
}
/*
//=============================================================================
//
//=============================================================================
void make_node_vector(Mesh* mesh)
{
        int ni;
        int nei;
        TNode* nod;
        ElementIter ele;
        double cnt;

        for( ni = 0; ni < mesh->n_nodes; ni++ ) {
                nod = mesh->node + ni;
                nod->vector[ 0 ] = 0.0;
                nod->vector[ 1 ] = 0.0;
                cnt = 0.0;
                for( nei = 0; nei < nod->n_elements(); nei++ ) {
                        ele = nod->element[ nei ];
                        nod->vector[ 0 ] += ele->vector[ 0 ];
                        nod->vector[ 1 ] += ele->vector[ 1 ];
                        cnt += 1.0;
                }
                nod->vector[ 0 ] /= cnt;
                nod->vector[ 1 ] /= cnt;

        }
}
 */
//=============================================================================
// FILL TH "FLUX" FIELD FOR ALL VV NEIGHBOURS IN THE MESH
//=============================================================================

void make_neighbour_flux(struct Problem *problem, Mesh* mesh) {
    struct Neighbour *ngh;

    FOR_NEIGHBOURS(ngh) {
        if (ngh->type != VV_2E)
            continue;
        ngh->flux = ngh->sigma * ngh->geom_factor * (ngh->element[1]->scalar - ngh->element[0]->scalar);
        continue;
    }
}
//=============================================================================
// FILL THE "PSCALAR" FIELD FOR ALL ELEMENTS IN THE MESH
//=============================================================================

void make_previous_scalar(struct Problem *problem, Mesh* mesh) {
    ElementIter ele;

    FOR_ELEMENTS(ele) {
        ele->pscalar = ele->scalar;
    }
}
//=============================================================================
//
//=============================================================================

void calc_element_balance(Mesh* mesh) {
    int li;
    ElementIter ele;

    xprintf(Msg, "Calculating balances... ")/*orig verb 2*/;

    FOR_ELEMENTS(ele) {
        ele->balance = 0.0;
        for (li = 0; li < ele->n_sides; li++)
            ele->balance += ele->side[ li ]->flux;
    }
    xprintf(Msg, "O.K.\n")/*orig verb 2*/;
}
//=============================================================================
//
//=============================================================================

void calc_external_balance(Mesh* mesh) {
    int bi;
    struct Boundary *bcd;
    double bal;

    xprintf(Msg, "Checking the external balance... ")/*orig verb 2*/;
    bal = 0.0;

    FOR_BOUNDARIES(bcd) {
        bal += bcd->side->flux;
    }
    if (fabs(bal) < ZERO)
        xprintf(Msg, "O.K.\n")/*orig verb 2*/;
    else
        xprintf(Msg, "Waters not balanced, error = %lg\n", bal)/*orig verb 2*/;
}
//=============================================================================
//
//=============================================================================

void water_balance(Mesh* mesh, MaterialDatabase* mat_base) {
    F_ENTRY;

    double bal;
    int c_water;
    struct Boundary *bcd;

    xprintf(Msg, "Calculating balance of water by types...\n")/*orig verb 2*/;
    FOR_BOUNDARIES(bcd)
    bcd->aux = 0;

    FOR_BOUNDARIES(bcd) {
        if (bcd->aux == 1)
            continue;
        c_water = bcd->group;
        bal = calc_water_balance(mesh, c_water);
        xprintf(Msg, "Water #%d:\t%g\n", c_water, bal);
    }


    xprintf(Msg, "Calculating sources of water by material types...\n");

    std::vector<double> source(mat_base->size(), 0); // initialize by zero
    MaterialDatabase::Iter mat;

    FOR_ELEMENTS(elm) {
        mat = mat_base->find_id(elm->mid);
        if (elm->source != NULL) {
            source[mat_base->index(mat)] += elm->volume * elm->source->density;
        }

    }

    FOR_MATERIALS_IT(*mat_base, mat) {
        xprintf(Msg, "Material #%d:\t% g\n", mat_base->get_id(mat), source[mat_base->index(mat)]);
    }


}
//=============================================================================
//
//=============================================================================

double calc_water_balance(Mesh* mesh, int c_water) {
    double rc;
    struct Boundary *bcd;

    rc = 0.0;

    FOR_BOUNDARIES(bcd) {
        if (bcd->aux == 1 || bcd->group != c_water)
            continue;
        rc += bcd->side->flux;
        bcd->aux = 1;
    }
    return rc;
}
