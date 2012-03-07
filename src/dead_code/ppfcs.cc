/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @ingroup io
 * @brief  Cross-sections computation
 *
 */

#include "transport/transport.h"

#include "system/system.hh"
#include "xio.h"
#include "system/math_fce.h"
//#include "problem.h"
#include "mesh/mesh.h"
#include "ppfcs.h"
#include "io/read_ini.h"
#include "materials.hh"


static void flow_cs(struct Transport *transport);
static void output_AGE(struct Transport *transport,double time);
static double particle_test(struct Transport *transport);
static void clear_tbc(struct Transport *transport);

static int create_flow_section(struct Transport *transport);
static void sfsec(struct Transport *transport);
static double node_pos(Node* node, double *eqn);
static void compute_cutplane(struct Transport *trannsport);
static double *cut_point(Node* node0, Node* node1, double eqn[4]);
static double *assign_node_coordinates(Node* node);
static struct ElementCut *new_ec(struct Transport *transport,ElementIter elm,struct ElementCut *prev_ec);
static void compute_flux(struct ElementCut *ec,struct FSection *fs);
static void output_FCS(struct Transport *transport);

//==============================================================================
//      Flow crosssection
//==============================================================================
void flow_cs(struct Transport *transport)
{
        struct ElementCut *ec;

        if (create_flow_section(transport)) {
                sfsec(transport);
                if(transport->fsec->elc != NULL){
                        compute_cutplane(transport);
                        FOR_ELEMENTCUT(ec)
                                compute_flux(ec,transport->fsec);
                        output_FCS(transport);
                }
        }
}
//==============================================================================
//      Create flow crosssection
//==============================================================================
int create_flow_section(struct Transport *transport)
{

    // This is ancient initialization from problem.c
    // TODO: Proper implementation of cross section
#if 0
    problem->ftrans_out       = get_b( "Output", "Write_ftrans_out", false );
    problem->cross_section    = get_b( "Output", "Cross_section", false );         //jh
    problem->cs_params        = get_s( "Output", "Cs_params", "0 0 0 0 0 0 0" );        //jh
//    problem->res_run          = get_b( "Output", "Cs_results_run", false );           //jh
//    problem->res_fin          = get_b( "Output", "Cs_results_final", false );
    problem->specify_elm_output =  get_b( "Output", "Specify_elm_type", false );   //jh temp
    problem->output_elm_type  = get_i( "Output", "Output_elm_type", 1 );        //jh temp
    problem->fsec_params       = get_s( "Output", "FCs_params", "0 0 0 0 0" );
//    problem->CF_params         = get_s( "Output", "ConfFlow_params", "0");
#endif


        double norm;
        int i;

        struct FSection *fsec;//=problem->fsec;

        transport->fsec =(struct FSection*)xmalloc(sizeof(struct FSection));
        fsec = transport->fsec;
        fsec->n_elm = 0;
        fsec->elc = NULL;

        fsec->fcs_params         = OptGetStr( "Output", "FCs_params", "0 0 0 0 0");

        if( sscanf( transport->fsec->fcs_params, "%lf %lf %lf %lf %d",
                &fsec->eqn[0],
                &fsec->eqn[1],
                &fsec->eqn[2],
                &fsec->eqn[3],
                &fsec->axis_output) == 5 )
            xprintf( Msg, "O.K.\n")/*orig verb 2*/;   //118
        else
            xprintf(UsrErr,"Invalid flow cross section params.\n");
        norm = vector_length(transport->fsec->eqn);

        if(norm > ZERO)
                for(i=0;i<4;i++)
                        fsec->eqn[i] /= norm;
        else{
                xfree(fsec);
                return false;
        }
        return true;
}
//==============================================================================
//      Select elements of flow crosssection
//==============================================================================
void sfsec(struct Transport *transport)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    ElementIter elm;
        NodeIter nod;
        struct FSection *fs;
        struct ElementCut *ec;
        int i,j;
        int b1,b2;

        fs = transport->fsec;
        ec = NULL;

        FOR_NODES( nod ){
                nod->faux = node_pos(nod,fs->eqn);
                nod->aux = SGN(nod->faux);
        }

        FOR_ELEMENTS(elm){
                b1 = b2 = 0;
                elm->aux = 0;
                FOR_ELEMENT_NODES(elm,i){ //element cut
                     if(elm->node[i]->aux == -1){
                        b1 = 1;
                        elm->aux = -1;
                    }
                     if(elm->node[i]->aux == 1)
                        b2 = 1;
                }
                if((b1 * b2) != 0){
                        ec = new_ec(transport,elm,ec);
                        ec->type = 2;
                }
                else
                FOR_ELEMENT_SIDES(elm,i){     //side cut
                        b1 = 1;
                        FOR_SIDE_NODES(elm->side[i],j)
                                if(elm->side[i]->node[j]->aux != 0)
                                        b1 = 0;
                                if((b1 == 1) && (elm->aux == -1)){
                                        ec = new_ec(transport,elm,ec);
                                        ec->type = 1;
                                        ec->sid = i;
                                }
                }
        }
}
//==============================================================================
//      Create new element cut
//==============================================================================
struct ElementCut *new_ec(struct Transport *transport,ElementIter elm,struct ElementCut *prev_ec)
{
        struct ElementCut *ec;
        ec = (struct ElementCut*)xmalloc(sizeof(struct ElementCut));

        ec->n_point = 0;
        ec->point[0] = NULL;
        ec->point[1] = NULL;
        ec->point[2] = NULL;
        ec->point[3] = NULL;

        if(transport->fsec->elc == NULL)
                transport->fsec->elc = ec;
        else
                prev_ec->next = ec;
        ec->prev = prev_ec;
        ec->next = NULL;
        transport->fsec->n_elm++;
        ec->element = elm;
        return ec;
}
//==============================================================================
//      Node position
//==============================================================================
double node_pos(Node* node, double* eqn)
{
        double val;

        val = eqn[0] * node->getX() +
              eqn[1] * node->getY() +
              eqn[2] * node->getZ() +
              eqn[3];

        return val;
}
//==============================================================================
//      Compute breakplane
//==============================================================================
void compute_cutplane(struct Transport *transport)
{
        ElementIter elm;
        struct ElementCut *ec;
        int i,j,n;

        FOR_ELEMENTCUT(ec){
                elm = ec->element;
                n = 0;
                FOR_ELEMENT_NODES(elm,i)
                        if(elm->node[i]->aux == 0)
                                ec->point[ec->n_point++] = assign_node_coordinates(elm->node[i]);
                if(ec->type != 1)
                        FOR_ELEMENT_NODES(elm,i){
                                n++;
                                for(j=n;j<elm->n_nodes;j++)
                                        if((elm->node[i]->aux != 0) && (elm->node[j]->aux != 0)){
                                                ec->point[ec->n_point] =
                                                cut_point(elm->node[i],elm->node[j],transport->fsec->eqn);
                                                if(ec->point[ec->n_point] != NULL)
                                                        ec->n_point++;
                                        }
                        }
        }

     /*   FOR_ELEMENTPLANE(ep){
                printf("\n\[id:%d\tdim:%d\tn_point:%d\ttype:%d\]",ep->element->id,ep->element->dim,ep->n_point,ep->type);
                for(j=0;j<ep->n_point;j++){
                        printf("\n");
                        for(k=0;k<3;k++)
                                printf("\t%f",ep->point[j][k]);
                }
                printf("\nBalance:%f\n",ep->element->balance);
        }
        printf("\nn_elm:%d",problem->fsec->n_elm);
        printf("\n%s",problem->fsec_params);
        getchar();

*/
}
//==============================================================================
//      Compute break point
//==============================================================================
double *cut_point(Node* node0, Node* node1, double eqn[4])
{
        double v[3],t,d,*point;

                point = NULL;

                v[0] = node0->getX() - node1->getX();
                v[1] = node0->getY() - node1->getY();
                v[2] = node0->getZ() - node1->getZ();
                t = node1->faux;
                d = scalar_product(eqn,v);
                if (fabs(d) < ZERO)
                        return point;
                else
                        t /= -d;

                if((t <= 1) && (t > ZERO)){
                        point = (double*)xmalloc(3*sizeof(double));
                        point[0] = node1->getX() + t * v[0];
                        point[1] = node1->getY() + t * v[1];
                        point[2] = node1->getZ() + t * v[2];
                }
                return point;
}
//==============================================================================
//      Assign node coordinates
//==============================================================================
double *assign_node_coordinates(Node* node)
{
        double *point;

                point = (double*)xmalloc(3*sizeof(double));
                point[0] = node->getX();
                point[1] = node->getY();
                point[2] = node->getZ();

                return point;
}
//==============================================================================
//      Compute flux
//==============================================================================
void compute_flux(struct ElementCut *ec,struct FSection *fs)
{
        ElementIter elm;
        double vector[3],u[3],v[3],en[3],flux,l;
        int i,sg;


        for(i=0;i<3;i++)
                vector[i] = 0.0;
        elm = ec->element;

        if(ec->type != 1)
                switch(elm->dim){
                        case 1:
                           /*
                        l = line_length(elm->node[1]->getX(),elm->node[1]->getY(),elm->node[1]->getZ(),elm->node[0]->getX(),elm->node[0]->getY(),elm->node[0]->getZ());

                                if(elm->side[0]->aux == 1){
                                        ld = line_length(ep->point[0][0],ep->point[0][1],ep->point[0][2],elm->node[0]->getX(),elm->node[0]->getY(),elm->node[0]->getZ());
                                        flux = elm->side[0]->flux + ld / l * (-elm->side[1]->flux - elm->side[0]->flux);
                                }
                                else{
                                        ld = line_length(ep->point[0][0],ep->point[0][1],ep->point[0][2],elm->node[1]->getX(),elm->node[1]->getY(),elm->node[1]->getZ());
                                        flux = elm->side[1]->flux + ld / l * (-elm->side[0]->flux - elm->side[1]->flux);
                                }
                                 */
                                flux = SGN(scalar_product(fs->eqn,elm->vector)) * elm->side[0]->metrics * elm->v_length;
                                ec->cutflux = flux;
                                break;
                        case 2:
                                u[ 0 ] = elm->node[ 1 ]->getX() - elm->node[ 0 ]->getX();
                                u[ 1 ] = elm->node[ 1 ]->getY() - elm->node[ 0 ]->getY();
                                u[ 2 ] = elm->node[ 1 ]->getZ() - elm->node[ 0 ]->getZ();
                                v[ 0 ] = elm->node[ 2 ]->getX() - elm->node[ 0 ]->getX();
                                v[ 1 ] = elm->node[ 2 ]->getY() - elm->node[ 0 ]->getY();
                                v[ 2 ] = elm->node[ 2 ]->getZ() - elm->node[ 0 ]->getZ();
                                vector_product( u, v, en );
                                normalize_vector( en );
                                vector_difference(ec->point[0],ec->point[1],v);
                                vector_product( v, en, u );
                                normalize_vector(u);
                                l = vector_length(v);
                                sg = SGN(scalar_product(elm->vector,fs->eqn));
                                ec->cutflux = fabs(scalar_product(u,elm->vector)) * elm->material->size * l * sg;
                                //ec->cutflux *= elm->material->size * l * sg;
                                break;
                        case 3:
                                break;
        }
        else
                ec->cutflux = elm->side[ec->sid]->flux;
}

//==============================================================================
//      OUTPUT FLOW CROSSSECTION
//==============================================================================
void output_FCS(struct Transport *transport)
{
        FILE *out;
        struct ElementCut *ec;
        double flux = 0;

	char dbl_fmt[ 16 ],file[LINE_SIZE];

 	sprintf( dbl_fmt, "%%.%dg", ConstantDB::getInstance()->getInt("Out_digit"));
        sprintf(file,"%s.fcs", IONameHandler::get_instance()->get_output_file_name(OptGetFileName("Output", "Output_file", NULL)) );
        out = xfopen(file,"wt");

        FOR_ELEMENTCUT(ec){
                flux += ec->cutflux;
                switch(ec->element->dim){
                 case 1:
          //              xfprintf(out,"%d\t",ec->element->id);
                        xfprintf(out,dbl_fmt,ec->point[0][transport->fsec->axis_output]);
                        xfprintf(out,"\t");
                        xfprintf(out,dbl_fmt,ec->cutflux);
                        break;
                 case 2:
                        if(ec->point[0][transport->fsec->axis_output] > ec->point[1][transport->fsec->axis_output]){
                          xfprintf(out,dbl_fmt,ec->point[0][transport->fsec->axis_output]);
                          xfprintf(out,":");
                          xfprintf(out,dbl_fmt,ec->point[1][transport->fsec->axis_output]);
                        }
                        else{
                          xfprintf(out,dbl_fmt,ec->point[1][transport->fsec->axis_output]);
                          xfprintf(out,":");
                          xfprintf(out,dbl_fmt,ec->point[0][transport->fsec->axis_output]);
                        }
                        xfprintf(out,"\t");
                        xfprintf(out,dbl_fmt,ec->cutflux);
                        break;
                }
                xfprintf(out,"\n");
        }

        xfprintf(out,"Total flux of cut : ");
        xfprintf(out,dbl_fmt,flux);
      //  printf("\nTotal flux of cut: %f\n",flux);
        xfprintf(out,"\n");
        xfclose( out );

}
//==============================================================================

// DECOVALEX

//==============================================================================
/*
//==============================================================================
//      OUTPUT AGE
//==============================================================================
void output_AGE(struct Transport *transport,double time)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    FILE *out;
    struct ElementCut *ec;
    double Np = 0.0;
    double pc = 0.0;
    double flux;
    char dbl_fmt[ 16 ],file[LINE_SIZE];
    int i,id;

    sprintf( dbl_fmt, "%%.%dg", ConstantDB::getInstance()->getInt("Out_digit"));
    sprintf(file,"%s.age", OptGetFileName("Output", "Output_file", NULL) );
    out = xfopen(file,(time==0.0) ? "wt":"at");

    // TODO: remove epos_id and find_id
        for(i=0;i<mesh->n_elements();i++)
                pc += transport->conc[MOBILE][0][i] * mesh->element.find_id(mesh->epos_id[i])->volume;

        FOR_ELEMENTCUT(ec){
        	id = id2pos(mesh,ELEMENT_FULL_ITER(ec->element).id(),mesh->epos_id,ELM);
        	flux = ec->cutflux;
        	Np +=  transport->conc[MOBILE][0][id] * flux * transport->time_step;
			}


        if((Np > 0.0) || (time == 0.0)){
         xfprintf(out,dbl_fmt,time);
         xfprintf(out,"\t");
         xfprintf(out,dbl_fmt,Np);
         xfprintf(out,"\t");
         xfprintf(out,dbl_fmt,pc);
         xfprintf(out,"\n");
         }

        xfclose( out );
       // getchar();
} */
//==============================================================================
//      PARTICLE TEST
//==============================================================================
/*
double particle_test(struct Transport *transport)
{
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    int i;
	double pc = 0.0;

	for(i=0;i < mesh->n_elements();i++){
		pc += transport->conc[MOBILE][0][i] * mesh->element.find_id(mesh->epos_id[i])->volume;
	}
	return pc;
} */
//==============================================================================
//      CLEAR TRANSPORT BC
//==============================================================================
void clear_tbc(struct Transport *transport)
{
        Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

	int start,stop,i;
	start = mesh->n_elements();
	stop = mesh->n_elements() + mesh->n_boundaries(); //-1
	for(i=start;i<stop;i++){
		transport->conc[MOBILE][0][i] = 0.0;
		transport->pconc[MOBILE][0][i] = 0.0;
	}
}
//------------------------------------------------------------------------------
