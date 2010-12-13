/*
 * convert.cc
 *
 *  Created on: 1.6.2010
 *      Author: jiri
 */

#include "constantdb.h"
#include "mesh/ini_constants_mesh.hh"
#include "transport_bcd.h"
#include "transport.h"

#include "system.hh"
#include "xio.h"
#include "output.h"
#include "math_fce.h"
#include "mesh.h"
#include "boundaries.h"
#include "output.h"
#include "convert.h"

#include "problem.h"
#include "sources.h"
//#include "materials.hh"
#include "concentrations.h"
#include "boundaries.h"
#include "postprocess.h"
#include "neighbours.h"

static void output_convert_to_pos_source(struct Problem *problem);
static void output_convert_to_pos_bcd(struct Problem *problem);
static void output_convert_to_pos_material(struct Problem *problem);
static void output_convert_to_pos_concentration(struct Problem *problem);
static void output_convert_to_pos_transport_bcd(struct Problem *problem);
/*
static void output_transport_convert(struct Problem *problem);
static void write_transport_ascii_data(FILE *out, struct Problem *problem, struct TNode **nodes, struct TElement **elements, int time_steps, int ph);
static void write_transport_binary_data(FILE *out, struct Problem *problem, struct TNode **nodes, struct TElement **elements, int time_steps, int ph);
*/

//=============================================================================
// OUTPUT ROUTINE FOR CONVERTING TO POS
//=============================================================================
void output_convert_to_pos(struct Problem *problem)
{
    F_ENTRY;

    ASSERT(!( problem == NULL ),"NULL as argument of function output_convert_to_pos()\n");
    if( OptGetBool("Output", "Write_output_file", "no") == false )
        return;
    if (ConstantDB::getInstance()->getChar("Sources_fname") != NULL)
            output_convert_to_pos_source(problem);
    output_convert_to_pos_bcd(problem);
    output_convert_to_pos_material(problem);
    if (OptGetBool("Transport", "Transport_on", "no") == true)
    {
            if (ConstantDB::getInstance()->getChar("Concentration_fname") != NULL)
                    output_convert_to_pos_concentration(problem);
            if (ConstantDB::getInstance()->getChar("Transport_bcd_fname") != NULL)
                    output_convert_to_pos_transport_bcd(problem);
    }
}
//=============================================================================
// OUTPUT ROUTINE FOR CONVERTING SOURCES TO POS
//=============================================================================
void output_convert_to_pos_source(struct Problem *problem)
{
  Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

  FILE *out;
  int li;
  char* filename = (char*)ConstantDB::getInstance()->getChar("Sources_fname");
  char dbl_fmt[ 16 ];
  ElementIter elm;
  Node *nod;

  F_ENTRY;

  ASSERT(!( problem == NULL ),"NULL as argument of function output_convert_to_pos_source()\n");
  sprintf( dbl_fmt, "%%.%dg ", ConstantDB::getInstance()->getInt("Out_digit"));
  strcat(filename,".pos");
  out = xfopen( filename, "wt" );
  xfprintf( out, "View \"%s - sources\" {\n", OptGetStr("Global", "Description", "No description.") );
  FOR_ELEMENTS(elm)
  {
    switch( elm->type ) {
      case LINE:
        xfprintf( out, "SL (" );
      break;
      case TRIANGLE:
        xfprintf( out, "ST (" );
      break;
      case TETRAHEDRON:
        xfprintf( out, "SS (" );
      break;
    }
    for( li = 0; li < elm->n_nodes; li++ ) {
      nod = elm->node[ li ];
      xfprintf( out, dbl_fmt, nod->getX() );
      xfprintf( out, ", " );
      xfprintf( out, dbl_fmt, nod->getY() );
      xfprintf( out, ", " );
      xfprintf( out, dbl_fmt, nod->getZ() );
      if( li == elm->n_nodes - 1 )
        xfprintf( out, ") {" );
      else
        xfprintf( out, ", " );
    }
    for( li = 0; li < elm->n_nodes; li++ ) {
      if (elm->source != NULL)
        xfprintf( out, dbl_fmt, elm->source->density);
      else
        xfprintf( out, dbl_fmt, 0.0);
      if( li == elm->n_nodes - 1 )
        xfprintf( out, "};\n" );
      else
        xfprintf( out, ", " );
    }
  }
  xfprintf( out, "};\n" );
  xfclose(out);
}
//=============================================================================
// OUTPUT ROUTINE FOR CONVERTING BCDs TO POS
//=============================================================================
void output_convert_to_pos_bcd(struct Problem *problem)
{
  Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

  const int bcd_types[] = {DIRICHLET,NEUMANN,NEWTON,NEWTON};
  const char *bcd_names[] = {"Dirichlet","Neumann","Newton P","Newton Sigma"};
  FILE *out;
  int li,j,test;
  char *filename = OptGetStr( "Input", "Boundary", "\\" );
  char dbl_fmt[ 16 ];
  ElementIter elm;
  Node* nod;
  struct Boundary *bcd;

  F_ENTRY;

  ASSERT(!( problem == NULL ),"NULL as argument of function output_convert_to_pos_bcd()\n");
  sprintf( dbl_fmt, "%%.%dg ", ConstantDB::getInstance()->getInt("Out_digit"));
  strcat(filename,".pos");
  out = xfopen( filename, "wt" );
  xfprintf( out, "View \"%s - mesh\" {\n", OptGetStr("Global", "Description", "No description.") );
  FOR_ELEMENTS(elm)
  {
    switch( elm->type ) {
      case LINE:
        xfprintf( out, "SL (" );
      break;
      case TRIANGLE:
        xfprintf( out, "ST (" );
      break;
      case TETRAHEDRON:
        xfprintf( out, "SS (" );
      break;
    }
    for( li = 0; li < elm->n_nodes; li++ ) {
      nod = elm->node[ li ];
      xfprintf( out, dbl_fmt, nod->getX() );
      xfprintf( out, ", " );
      xfprintf( out, dbl_fmt, nod->getY() );
      xfprintf( out, ", " );
      xfprintf( out, dbl_fmt, nod->getZ() );
      if( li == elm->n_nodes - 1 )
        xfprintf( out, ") {" );
      else
        xfprintf( out, ", " );
    }
    for( li = 0; li < elm->n_nodes; li++ ) {
      xfprintf( out, dbl_fmt, 0.0);
      if( li == elm->n_nodes - 1 )
        xfprintf( out, "};\n" );
      else
        xfprintf( out, ", " );
    }
  }
  xfprintf( out, "};\n" );
  for (j = 0;j < 4; j++)
  {
    test = false;
    FOR_BOUNDARIES(bcd)
      if (bcd->type == bcd_types[ j ])
      {
        test = true;
        break;
      }
    if (test == false)
      continue;
    xfprintf( out, "View \"%s - bcd %s\" {\n", OptGetStr("Global", "Description", "No description."), bcd_names[j] );
    FOR_BOUNDARIES(bcd)
    {
      if (bcd->type != bcd_types[ j ]){
        continue;
      }
      switch( bcd->side->shape ) {
        case xPOINT:
            xfprintf( out, "SP (" );
        break;
        case LINE:
            xfprintf( out, "SL (" );
        break;
        case TRIANGLE:
            xfprintf( out, "ST (" );
        break;
      }
      for( li = 0; li < bcd->side->n_nodes; li++ ) {
        nod = bcd->side->node[ li ];
        xfprintf( out, dbl_fmt, nod->getX() );
        xfprintf( out, ", " );
        xfprintf( out, dbl_fmt, nod->getY() );
        xfprintf( out, ", " );
        xfprintf( out, dbl_fmt, nod->getZ() );
        if( li == bcd->side->n_nodes - 1 )
          xfprintf( out, ") {" );
        else
          xfprintf( out, ", " );
      }
      for( li = 0; li < bcd->side->n_nodes; li++ ) {
        switch (bcd->type)
        {
          case DIRICHLET: xfprintf( out, dbl_fmt, bcd->scalar);break;
          case NEUMANN: xfprintf( out, dbl_fmt, bcd->flux);break;
          case NEWTON:
            switch (j)
            {
              case 2: xfprintf( out, dbl_fmt, bcd->scalar);break;
              case 3: xfprintf( out, dbl_fmt, bcd->sigma);break;
            }
          break;
        }
        if( li == bcd->side->n_nodes - 1 )
          xfprintf( out, "};\n" );
        else
          xfprintf( out, ", " );
      }
    }
    xfprintf( out, "};\n" );
  }
  xfclose(out);
}
//=============================================================================
// OUTPUT ROUTINE FOR CONVERTING MATERIALS TO POS
//=============================================================================
void output_convert_to_pos_material(struct Problem *problem)
{
  Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

  FILE* out;
  int li;
  char* filename = (char*)ConstantDB::getInstance()->getChar("Material_fname");
  char dbl_fmt[ 16 ];
  char int_fmt[ 16 ];
  ElementIter elm;
  Node* nod;

  F_ENTRY;

  sprintf( dbl_fmt, "%%.%dg ", ConstantDB::getInstance()->getInt("Out_digit"));
 // sprintf( int_fmt, "%%d ", problem->out_digit );
  sprintf( int_fmt, "%d ", ConstantDB::getInstance()->getInt("Out_digit"));
  strcat( filename, ".pos" );
  out = xfopen( filename, "wt" );
  xfprintf( out, "View \"%s - materials\" {\n", OptGetStr("Global", "Description", "No description.") );
  FOR_ELEMENTS(elm)
  {
    switch( elm->type ) {
      case LINE:
        xfprintf( out, "SL (" );
      break;
      case TRIANGLE:
        xfprintf( out, "ST (" );
      break;
      case TETRAHEDRON:
        xfprintf( out, "SS (" );
      break;
    }
    for( li = 0; li < elm->n_nodes; li++ ) {
      nod = elm->node[ li ];
      xfprintf( out, dbl_fmt, nod->getX() );
      xfprintf( out, ", " );
      xfprintf( out, dbl_fmt, nod->getY() );
      xfprintf( out, ", " );
      xfprintf( out, dbl_fmt, nod->getZ() );
      if( li == elm->n_nodes - 1 )
        xfprintf( out, ") {" );
      else
        xfprintf( out, ", " );
    }
    for( li = 0; li < elm->n_nodes; li++ ) {
      xfprintf( out, int_fmt, elm->mid );
      if( li == elm->n_nodes - 1 )
        xfprintf( out, "};\n" );
      else
        xfprintf( out, ", " );
    }
  }
  xfprintf( out, "};\n" );
  xfclose(out);
}
//=============================================================================
// OUTPUT ROUTINE FOR CONVERTING CONCENTRATIONS TO POS
//=============================================================================
void output_convert_to_pos_concentration(struct Problem *problem)
{
  Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

  FILE *out;
  int li;
  char* filename = (char*)ConstantDB::getInstance()->getChar("Concentration_fname");
  char dbl_fmt[ 16 ];
  ElementIter elm;
  Node* nod;

  F_ENTRY;

  sprintf( dbl_fmt, "%%.%dg ", ConstantDB::getInstance()->getInt("Out_digit"));
  strcat(filename,".pos");
  out = xfopen( filename, "wt" );
  xfprintf( out, "View \"%s - Concentrations\" {\n", OptGetStr("Global", "Description", "No description.") );
  FOR_ELEMENTS(elm)
  {
    switch( elm->type ) {
      case LINE:
        xfprintf( out, "SL (" );
      break;
      case TRIANGLE:
        xfprintf( out, "ST (" );
      break;
      case TETRAHEDRON:
        xfprintf( out, "SS (" );
      break;
    }
    for( li = 0; li < elm->n_nodes; li++ ) {
      nod = elm->node[ li ];
      xfprintf( out, dbl_fmt, nod->getX() );
      xfprintf( out, ", " );
      xfprintf( out, dbl_fmt, nod->getY() );
      xfprintf( out, ", " );
      xfprintf( out, dbl_fmt, nod->getZ() );
      if( li == elm->n_nodes - 1 )
        xfprintf( out, ") {" );
      else
        xfprintf( out, ", " );
    }
    for( li = 0; li < elm->n_nodes; li++ ) {
      if (elm->start_conc != NULL)
        xfprintf( out, dbl_fmt, elm->start_conc->conc);
      else
        xfprintf( out, dbl_fmt, 0.0);
      if( li == elm->n_nodes - 1 )
        xfprintf( out, "};\n" );
      else
        xfprintf( out, ", " );
    }
  }
  xfprintf( out, "};\n" );
  xfclose(out);
}
//=============================================================================
// OUTPUT ROUTINE FOR CONVERTING TRANSPORT_BCDs TO POS
//=============================================================================
void output_convert_to_pos_transport_bcd(struct Problem *problem)
{
  Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

  FILE *out;
  int li;
  char *filename = (char*)ConstantDB::getInstance()->getChar("Transport_bcd_fname");
  char dbl_fmt[ 16 ];
  ElementIter elm;
  Node* nod;
  struct Boundary *bcd;

  F_ENTRY;

  sprintf( dbl_fmt, "%%.%dg ", ConstantDB::getInstance()->getInt("Out_digit"));
  strcat(filename,".pos");
  out = xfopen( filename, "wt" );
  xfprintf( out, "View \"%s - Mesh\" {\n", OptGetStr("Global", "Description", "No description.") );
  FOR_ELEMENTS(elm)
  {
    switch( elm->type ) {
      case LINE:
        xfprintf( out, "SL (" );
      break;
      case TRIANGLE:
        xfprintf( out, "ST (" );
      break;
      case TETRAHEDRON:
        xfprintf( out, "SS (" );
      break;
    }
    for( li = 0; li < elm->n_nodes; li++ ) {
      nod = elm->node[ li ];
      xfprintf( out, dbl_fmt, nod->getX() );
      xfprintf( out, ", " );
      xfprintf( out, dbl_fmt, nod->getY() );
      xfprintf( out, ", " );
      xfprintf( out, dbl_fmt, nod->getZ() );
      if( li == elm->n_nodes - 1 )
        xfprintf( out, ") {" );
      else
        xfprintf( out, ", " );
    }
    for( li = 0; li < elm->n_nodes; li++ ) {
      xfprintf( out, dbl_fmt, 0.0);
      if( li == elm->n_nodes - 1 )
        xfprintf( out, "};\n" );
      else
        xfprintf( out, ", " );
    }
  }
  xfprintf( out, "};\n" );
  xfprintf( out, "View \"%s - Transport Boundary Conditions\" {\n", OptGetStr("Global", "Description", "No description.") );
  FOR_BOUNDARIES(bcd)
  {
    switch( bcd->side->shape ) {
      case xPOINT:
        xfprintf( out, "SP (" );
      break;
      case LINE:
        xfprintf( out, "SL (" );
      break;
      case TRIANGLE:
        xfprintf( out, "ST (" );
      break;
    }
    for( li = 0; li < bcd->side->n_nodes; li++ ) {
      nod = bcd->side->node[ li ];
      xfprintf( out, dbl_fmt, nod->getX() );
      xfprintf( out, ", " );
      xfprintf( out, dbl_fmt, nod->getY() );
      xfprintf( out, ", " );
      xfprintf( out, dbl_fmt, nod->getZ() );
      if( li == bcd->side->n_nodes - 1 )
        xfprintf( out, ") {" );
      else
        xfprintf( out, ", " );
    }
    for( li = 0; li < bcd->side->n_nodes; li++ ) {
      xfprintf( out, dbl_fmt, bcd->transport_bcd->conc);
      if( li == bcd->side->n_nodes - 1 )
        xfprintf( out, "};\n" );
      else
        xfprintf( out, ", " );
    }
  }
  xfprintf( out, "};\n" );
  xfclose(out);
}

// folowing function seems to be completly WRONG by desing and implementation
// - no need to recreate whole Element and Node structures HERE
// - wrong allocation
//
#if 0

//==============================================================================
//      OUTPUT FLOW CONVERT
//==============================================================================
void output_convert(struct Problem *problem){
    /*    struct TNode{
                double *scalar;
        };
        struct TElement{
                double *scalar;
                double **vector;
        };  */
        Mesh* mesh;
        ElementIter ele;
        Node* nod;
        struct Node *nodes;
        struct TElement *elements;
        FILE *out,*in;
    char dbl_fmt[ 16 ];
        char filenamei[255];
        char filenameo[255];
        char line[LINE_SIZE];
        int li,time_steps,j,id,k,ci;
        double value;
        bool found;

        F_ENTRY;

        mesh = problem->mesh;
      //  nodes = (struct TNode**)xmalloc(sizeof(struct TNode*));
      //  elements = (struct TElement**)xmalloc(sizeof(struct TElement*));
    sprintf( dbl_fmt, "%%.%dg ", problem->out_digit );
        sprintf( filenameo,"%s",problem->out_fname);
        sprintf( filenamei,"%s.tmp",problem->out_fname);
        out = xfopen( filenameo, "wt" );
        in = xfopen( filenamei, "rt" );

        found=skip_to(in, "VALUES");
        ASSERT(found, "Can not find section: VALUES\n");

    // find number of time steps
        time_steps = 0;
        while(skip_to(in, "Time"))
            time_steps++;

        nodes = (struct Node *)xmalloc((mesh->max_nod_id + 1)* sizeof(struct Node));
        elements = (struct TElement *)xmalloc((mesh->max_elm_id + 1)* sizeof(struct TElement));

    // allocate and initialize tmp nodes scalar
        for (li = 0; li <= mesh->max_nod_id; li++){
                nodes[li].scalar = (double *)xmalloc(time_steps* sizeof(double));
                for (j = 0; j < time_steps; j++)
                        nodes[li].scalar[j] = 0;
        }


    // allocate and initialize tmp element scalar & vector
        for (li = 0; li <= mesh->max_elm_id; li++){
                elements[li].scalar = (double *)xmalloc( time_steps * sizeof(double));
                elements[li].vector = (double **)xmalloc(3 * sizeof(double*));
        for( ci = 0; ci < 3; ci++ )
            elements[li].vector[ci] = (double *)xmalloc(time_steps * sizeof(double));
                for (j = 0; j < time_steps; j++)
                {
                        elements[li].scalar[j] = 0;
                for( ci = 0; ci < 3; ci++ )
                                elements[li].vector[ci][j] = 0;
                }
        }

    // now read values for nodes
        xfclose(in);
        in = xfopen( filenamei, "rt" );
        found=skip_to(in, "VALUES");
        ASSERT(found, "Can not find section: VALUES\n");

        j = 0;
        while(skip_to(in, "Nodes")){
                for (li = 0; li < mesh->n_nodes; li++){
                        xfgets( line, LINE_SIZE - 2, in );
                        id = atoi( xstrtok( line) );
                        value = atof ( xstrtok( NULL) );
                        nodes[id].scalar[j] = value;
                }
        j++;
        }

    // end of TMP file reading
    // read values for elements
        xfclose(in);
        in = xfopen( filenamei, "rt" );
        found=skip_to(in, "VALUES");
        ASSERT(found, "Can not find section: VALUES\n");

        j = 0;
        while(skip_to(in, "Elements")){
                for (li = 0; li < mesh->n_elements(); li++){
                        xfgets( line, LINE_SIZE - 2, in );
                        id = atoi( xstrtok( line) );
                        value = atof ( xstrtok( NULL) );
                        elements[id].scalar[j] = value;
                        for( ci = 0; ci < 3; ci++ ) {
                                value = atof ( xstrtok( NULL) );
                                elements[id].vector[ci][j] = value;
                        }
                }
                j++;
        }




        write_ascii_header(problem,out);

        for( k = 0; k < 4; k++){  // p node, p element, pz element, u element
        switch(k){
                  case 0:
                        xfprintf( out, "View \"%s - p\" {\n", problem->description );
                        break;
                  case 1:
                        xfprintf( out, "View \"%s - pc\" {\n", problem->description );
                        break;
                  case 2:
                        xfprintf( out, "View \"%s - pz\" {\n", problem->description );
                        break;
                  case 3:
                        xfprintf( out, "View \"%s - u\" {\n", problem->description );
                        break; }


        if(k != 3)
        {
        FOR_ELEMENTS( ele ){
            switch( ele->type ) {
                case LINE:
                    xfprintf( out, "SL (" );
                    break;
                case TRIANGLE:
                    xfprintf( out, "ST (" );
                    break;
                case TETRAHEDRON:
                    xfprintf( out, "SS (" );
                    break;
                        }
            for( li = 0; li < ele->n_nodes; li++ ) {
                nod = ele->node[ li ];
                xfprintf( out, dbl_fmt, nod->x );
                xfprintf( out, ", " );
                xfprintf( out, dbl_fmt, nod->y );
                xfprintf( out, ", " );
                xfprintf( out, dbl_fmt, nod->z );
                if( li == ele->n_nodes - 1 )
                    xfprintf( out, ") {" );
                else
                    xfprintf( out, ", " );
            }

            for( li = 0; li < time_steps; li++ ) {
                        for(j = 0; j < ele->n_nodes; j++){
                switch(k){
                  case 0:
                        nod = ele->node[ j ];
                        xfprintf( out, dbl_fmt,nodes[nod->id].scalar[li]);
                        break;
                  case 1:
                        xfprintf( out, dbl_fmt,elements[ele->id].scalar[li]);
                        break;
                  case 2:
                        xfprintf( out, dbl_fmt,elements[ele->id].scalar[li] + ele->centre[2]);
                        break;
                  }
                        if( li == time_steps - 1 && j == ele->n_nodes - 1)
                        xfprintf( out, "};\n" );
                                else
                                xfprintf( out, ", " );
                        }
            }


        } // for elements
        } // end if k != 3
        else

        FOR_ELEMENTS( ele ){
                xfprintf( out, "VP (" );
        xfprintf( out, dbl_fmt, ele->centre[ 0 ] );
        xfprintf( out, ", " );
        xfprintf( out, dbl_fmt, ele->centre[ 1 ] );
        xfprintf( out, ", " );
        xfprintf( out, dbl_fmt, ele->centre[ 2 ] );
                xfprintf( out, ") {" );
                for( li = 0; li < time_steps; li++ ){
                xfprintf( out, dbl_fmt, elements[ele->id].vector[ 0 ][li] );
        xfprintf( out, ", " );
        xfprintf( out, dbl_fmt, elements[ele->id].vector[ 1 ][li] );
        xfprintf( out, ", " );
        xfprintf( out, dbl_fmt, elements[ele->id].vector[ 2 ][li] );
                if( li == time_steps - 1)
                        xfprintf( out, "};\n" );
                else
                        xfprintf( out, ", " );
                }
                }
        xfprintf( out, "};\n" );
        } //p node, p element, pz element, u element

        xfclose(in);
        xfclose(out);
        xfree(nodes);
        xfree(elements);

}
#endif

#if 0
//==============================================================================
//      WRITE TRANSPORT ASCII DATA TO THE POS FILE
//==============================================================================
static void write_transport_ascii_data(FILE *out, struct Problem *problem, struct Node **nodes, struct TElement **elements, int time_steps, int ph)
{
        int k,sbi,li,j,n_subst;
        Mesh* mesh;
        ElementIter ele;
        Node* nod;
        char dbl_fmt[ 16 ];
        double norm,vconc[3];

        n_subst = problem->transport->n_substances;
        mesh = problem->mesh;
        sprintf( dbl_fmt, "%%.%dg ", problem->out_digit );

        write_ascii_header(problem, out); // header in POS for user view

        for( k = 0; k < 2; k++){  // nodes and elements
    for( sbi = 0; sbi < n_subst; sbi++ ) {

                xfprintf( out, "View \"Concentration in %s of %s\" {\n",
                (k == 0) ? "node" : "element",problem->transport->substance_name[ sbi ] );


        FOR_ELEMENTS( ele ){
            switch( ele->type ) {
                case LINE:
                    xfprintf( out, "SL (" );
                    break;
                case TRIANGLE:
                    xfprintf( out, "ST (" );
                    break;
                case TETRAHEDRON:
                    xfprintf( out, "SS (" );
                    break;
                        }
                FOR_ELEMENT_NODES(ele,li){
                nod = ele->node[ li ];
                xfprintf( out, dbl_fmt, nod->x );
                xfprintf( out, ", " );
                xfprintf( out, dbl_fmt, nod->y );
                xfprintf( out, ", " );
                xfprintf( out, dbl_fmt, nod->z );
                if( li == ele->n_nodes - 1 )
                    xfprintf( out, ") {" );
                else
                    xfprintf( out, ", " );
            }

            for( li = 0; li < time_steps; li++ ) {
                        FOR_ELEMENT_NODES(ele,j){
                switch(k){
                  case 0:
                        nod = ele->node[ j ];
                        xfprintf( out, dbl_fmt,nodes[ph][nod->id].conc[sbi][li]);
                        break;
                  case 1:
                        xfprintf( out, dbl_fmt,elements[ph][ele->id].conc[sbi][li]);
                        break;
                  }
                        if( li == time_steps - 1 && j == ele->n_nodes - 1)
                        xfprintf( out, "};\n" );
                                else
                                xfprintf( out, ", " );
                        }
            }

        } // END FOR_ELEMENTS
        xfprintf( out, "};\n" );

       // CONC VECTOR IN MOBILE ZONE FOR STEADY SATURATED PROBEM
        if((ph == 0) & (k == 1) & (problem->type == STEADY_SATURATED)){
       xfprintf( out, "View \"Concentration vector of %s\" {\n",
            problem->transport->substance_name[ sbi ] );
         FOR_ELEMENTS( ele ){
                xfprintf( out, "VP (" );
        xfprintf( out, dbl_fmt, ele->centre[ 0 ] );
        xfprintf( out, ", " );
        xfprintf( out, dbl_fmt, ele->centre[ 1 ] );
        xfprintf( out, ", " );
        xfprintf( out, dbl_fmt, ele->centre[ 2 ] );
                xfprintf( out, ") {" );
                norm = vector_length(ele->vector);
                if(norm > 0){
                        vconc[0] = ele->vector[0]/norm; vconc[1] = ele->vector[1]/norm; vconc[2] = ele->vector[2]/norm;
                        }
                else
                        scale_vector(vconc,0);
                for( li = 0; li < time_steps; li++ ){
                xfprintf( out, dbl_fmt, elements[ph][ele->id].conc[sbi][li] * vconc[0]);
        xfprintf( out, ", " );
        xfprintf( out, dbl_fmt, elements[ph][ele->id].conc[sbi][li] * vconc[1]);
        xfprintf( out, ", " );
        xfprintf( out, dbl_fmt, elements[ph][ele->id].conc[sbi][li] * vconc[2]);
                if( li == time_steps - 1)
                        xfprintf( out, "};\n" );
                else
                        xfprintf( out, ", " );
                }
                }
        xfprintf( out, "};\n" );
        }
       // END CONC VECTOR
    } // end for substances
        } // end nodes and elements


}
//==============================================================================
//      WRITE TRANSPORT BINARY DATA TO THE POS FILE
//==============================================================================
static void write_transport_binary_data(FILE *out, struct Problem *problem, struct Node **nodes, struct TElement **elements, int time_steps, int ph)
{
        int k,sbi,i,li,j,n_subst;
        double ts;
        Mesh* mesh;
        ElementIter ele;
        int one = 1;
        double norm,vconc[3],vconct[3];

        n_subst = problem->transport->n_substances;
        mesh = problem->mesh;

        xfprintf(out, "$PostFormat\n");
        xfprintf(out, "%g %d %d\n", 1.4, 1, sizeof(double));
        xfprintf(out, "$EndPostFormat\n");

        for( k = 0; k < 2; k++){  // nodes and elements
            for( sbi = 0; sbi < n_subst; sbi++ ) {
                        xfprintf( out, "$View\nConcentration^in^%s^of^%s ",
            (k == 0) ? "node" : "element",problem->transport->substance_name[ sbi ] );
        //uprava
       // mesh->n_tetrahedras = 0;
        xfprintf(out, "%d 0 0 0 %d 0 0 %d 0 0 0 0 0 %d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \n",time_steps,mesh->n_lines,mesh->n_triangles,mesh->n_tetrahedras);
        xfwrite(&one, sizeof(int), 1, out);
        ts = 0;
        for(li=0;li<time_steps;li++){
                xfwrite(&ts, sizeof(double), 1, out);
                ts += 1;//problem->time_step;
                }
        for(i = 1;i < 4;i++)
                FOR_ELEMENTS( ele )
                        if(ele->dim == i){ // && ele->dim != 3){
                                write_elm_position_to_binary_output(ele,out);
                                for( j = 0; j < time_steps; j++ )
                                        FOR_ELEMENT_NODES(ele,li)
                                                xfwrite( (k==0) ?  &nodes[ph][ele->node[li]->id].conc[sbi][j] :
                                                &elements[ph][ele->id].conc[sbi][j],sizeof(double),1,out);
                        }
        xfprintf(out, "\n$EndView\n");
                }
        }
        //  Vector conc
        for( sbi = 0; sbi < n_subst; sbi++ ) {
                xfprintf( out, "$View\nConcentration^vector^of^%s ",problem->transport->substance_name[ sbi ] );
        xfprintf(out, "%d 0 %d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n",time_steps,mesh->n_elements());
        xfwrite(&one, sizeof(int), 1, out);
        ts = 0;
        for(li=0;li<time_steps;li++){
                xfwrite(&ts, sizeof(double), 1, out);
                ts += 1;//problem->time_step;
                }
        FOR_ELEMENTS( ele ){
                xfwrite(&ele->centre[0],sizeof(double),3,out);
                norm = vector_length(ele->vector);
                if(norm > 0){
                        vconc[0] = ele->vector[0]/norm; vconc[1] = ele->vector[1]/norm; vconc[2] = ele->vector[2]/norm;
                        }
                else
                        scale_vector(vconc,0);
                for( j = 0; j < time_steps; j++ ){
                        vconct[0] = vconc[0] * elements[ph][ele->id].conc[sbi][j];
                        vconct[1] = vconc[1] * elements[ph][ele->id].conc[sbi][j];
                        vconct[2] = vconc[2] * elements[ph][ele->id].conc[sbi][j];
                        xfwrite(&vconct[0],sizeof(double),3,out);
                }
         }
        xfprintf(out, "\n$EndView\n");
        }
}

//==============================================================================
//      OUTPUT TRANSPORT CONVERT
//==============================================================================
void output_transport_convert(struct Problem *problem){

        Mesh* mesh = problem->mesh;
        struct Node** nodes;
        struct TElement** elements;
        FILE **out,**in;
        char line[LINE_SIZE];
        int li,time_steps,j,id,ph;
        double value;
        int sbi;
        int n_subst = problem->transport->n_substances;
        bool found;


        nodes = (struct Node**)xmalloc((4)* sizeof(struct Node*));
        elements = (struct TElement**)xmalloc((4)* sizeof(struct TElement*));


        in = open_temp_files(problem->transport, "%s.tmp", "rt" );
        out = open_temp_files(problem->transport, "%s", (problem->pos_format_id == POS_ASCII) ? "wt" : "wb" );


    // find number of time steps
        found=skip_to(in[0], "VALUES");
        ASSERT(found , "Can not find section: VALUES\n");
        time_steps = 0;
        while(skip_to(in[0], "Time"))
                time_steps++;

        for(ph = 0; ph < 4; ph++){           // phase for cycle

        if (in[ph] != NULL){
                nodes[ph] = (struct Node *)xmalloc((mesh->max_nod_id + 1)* sizeof(struct Node));
                elements[ph] = (struct TElement *)xmalloc((mesh->max_elm_id + 1)* sizeof(struct TElement));
        }
        else
        continue;

    // allocate and initialize tmp nodes conc
        for (li = 0; li <= mesh->max_nod_id; li++){
                nodes[ph][li].conc = (double **)xmalloc(n_subst * sizeof(double*));
        for( sbi = 0; sbi < n_subst; sbi++ )
            nodes[ph][li].conc[sbi] = (double *)xmalloc(time_steps* sizeof(double));
        for( sbi = 0; sbi < n_subst; sbi++ )
                    for (j = 0; j < time_steps; j++)
                          nodes[ph][li].conc[sbi][j] = 0;
        }
    // allocate and initialize tmp element conc
        for (li = 0; li <= mesh->max_elm_id; li++){
                elements[ph][li].conc = (double **)xmalloc(n_subst * sizeof(double*));
        for( sbi = 0; sbi < n_subst; sbi++ )
            elements[ph][li].conc[sbi] = (double *)xmalloc(time_steps* sizeof(double));
        for( sbi = 0; sbi < n_subst; sbi++ )
                    for (j = 0; j < time_steps; j++)
                          elements[ph][li].conc[sbi][j] = 0;
        }

    // now read values for nodes
        xfclose(in[ph]);
        in = open_temp_files(problem->transport, "%s.tmp", "rt" );
        found=skip_to(in[ph], "VALUES");
        ASSERT(found, "Can not find section: VALUES\n");


        j = 0;
        while(skip_to(in[ph], "Nodes")){
            for (li = 0; li < mesh->n_nodes; li++){
                                xfgets( line, LINE_SIZE - 2, in[ph] );
                id = atoi( xstrtok( line) );
                for( sbi = 0; sbi < n_subst; sbi++ ) {
                    value = atof ( xstrtok( NULL) );
                    nodes[ph][id].conc[sbi][j] = value;
                }
             }
        j++;
        }


    // end of TMP file reading
    // read values for elements
        xfclose(in[ph]);
        in = open_temp_files(problem->transport, "%s.tmp", "rt" );
        found=skip_to(in[ph], "VALUES");
        ASSERT(found, "Can not find section: VALUES\n");

        j = 0;
        while(skip_to(in[ph], "Elements")){
            for (li = 0; li < mesh->n_elements(); li++){
                                xfgets( line, LINE_SIZE - 2, in[ph] );
                id = atoi( xstrtok( line) );
                for( sbi = 0; sbi < n_subst; sbi++ ) {
                    value = atof ( xstrtok( NULL) );
                    elements[ph][id].conc[sbi][j] = value;
                }
            }
                j++;
        }
        xfclose(in[ph]);

        xprintf( Msg, "Writing transport output files... ")/*orig verb 2*/;
                switch(problem->pos_format_id){
                        case POS_ASCII:
                                write_transport_ascii_data(out[ph],problem,nodes,elements,time_steps,ph);
                                break;
                        case POS_BIN:
                                write_transport_binary_data(out[ph],problem,nodes,elements,time_steps,ph);
                                break;
                }
                xfclose(out[ph]);
        }// end ph
        xprintf( Msg, "O.K.\n")/*orig verb 2*/;

        xfree(nodes);
        xfree(elements);
        xfree(in);
        xfree(out);

}
#endif

