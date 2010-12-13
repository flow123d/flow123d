/* 
 * File:   gmshmeshreader.cc
 * Author: dalibor
 * 
 * Created on October 3, 2010, 11:32 AM
 */

#include "msh_gmshreader.h"
#include "nodes.h"
#include "xio.h"

GmshMeshReader::GmshMeshReader() {
    xprintf(Msg, " - GmshMeshReader()\n");
}

GmshMeshReader::~GmshMeshReader() {
}

/**
 *  Read mesh from file
 */
void GmshMeshReader::read(const char* fileName, Mesh* mesh) {
    xprintf(Msg, " - GmshMeshReader->read(const char* fileName, Mesh* mesh)\n");

    ASSERT(!(fileName == NULL), "Argument fileName is NULL in method GmshMeshRedaer->read(const char*, Mesh*)\n");
    ASSERT(!(mesh == NULL), "Argument mesh is NULL in method GmshMeshRedaer->read(const char*, Mesh*)\n");

    FILE* file = xfopen(fileName, "rt");

    read_nodes(file, mesh);
    read_elements(file, mesh);

    xfclose(file);
}

/**
 * private method for reading of nodes
 */
void GmshMeshReader::read_nodes(FILE* in, Mesh* mesh) {
    xprintf(Msg, " - Reading nodes...");

    char line[ LINE_SIZE ];

    skip_to(in, "$Nodes");
    xfgets(line, (LINE_SIZE - 2), in);
    int numNodes = atoi(xstrtok(line));
    ASSERT(!(numNodes < 1), "Number of nodes < 1 in function read_node_list()\n");

    mesh->node_vector.reserve(numNodes);

    for (int i = 0; i < numNodes; ++i) {
        xfgets(line, LINE_SIZE - 2, in);

        int id = atoi(xstrtok(line));
        //TODO: co kdyz id <= 0 ???
        INPUT_CHECK(!(id < 0), "Id of node must be > 0\n");

        double x = atof(xstrtok(NULL));
        double y = atof(xstrtok(NULL));
        double z = atof(xstrtok(NULL));

        NodeFullIter node = mesh->node_vector.add_item(id);
        node->id = id;
        node->set(x,y,z);
    }

    //xprintf(MsgVerb, " %d nodes readed. ", nodeList->size());
    xprintf(Msg, " %d nodes readed. ", mesh->node_vector.size());
    xprintf(Msg, "O.K.\n");
}

/**
 * supported_element_type(int type)
 */
char GmshMeshReader::supported_element_type(int type) {
    switch (type) {
        case LINE:
        case TRIANGLE:
        case TETRAHEDRON:
            return true;
    }
    return false;
}

/**
 * PARSE ELEMENT LINE
 */
void GmshMeshReader::parse_element_line(ElementVector &ele_vec, char *line, Mesh* mesh) {
    int id, ti, ni;
    int type;
    int n_tags = NDEF;

    F_ENTRY;
    ASSERT(NONULL(line), "NULL as argument of function parse_element_line()\n");

    //get element ID
    id = atoi(xstrtok(line));
    INPUT_CHECK(id >= 0, "Id number of element must be >= 0\n");

    //DBGMSG("add el: %d", id);
    ElementFullIter ele(ele_vec.add_item(id));
    //get element type: supported:
    //	1 Line (2 nodes)
    //	2 Triangle (3 nodes)
    //	4 Tetrahedron (4 nodes)

    type = atoi(xstrtok(NULL));
    if (supported_element_type(type) == false)
        xprintf(UsrErr, "Element %d is of the unsupported type %d\n", id, type);

    ele->type = type;
    element_type_specific(ele);
    element_allocation_independent(ele);

    //get number of tags (at least 2)
    n_tags = atoi(xstrtok(NULL));
    INPUT_CHECK(!(n_tags < 2), "At least two element tags have to be defined. Elm %d\n", id);

    //get tags 1 and 2
    ele->mid = atoi(xstrtok(NULL));
    ele->rid = atoi(xstrtok(NULL));

    //get remaining tags
    if (n_tags > 2) {
        ele->pid = atoi(xstrtok(NULL)); // chop partition number in the new GMSH format

        //skip remaining tags
        for (ti = 3; ti < n_tags; ti++)
            xstrtok(NULL);
    }

    ele->node = (Node**) xmalloc(ele->n_nodes * sizeof (Node*));

    FOR_ELEMENT_NODES(ele, ni) {
        int nodeId = atoi(xstrtok(NULL));
        NodeIter node = mesh->node_vector.find_id( nodeId );

        ASSERT(NONULL(node), "Unknown node with label %d\n", nodeId);

        ele->node[ni] = node;
    }
}

/**
 * element_type_specific(ElementFullIter ele)
 */
void GmshMeshReader::element_type_specific(ElementFullIter ele) {
    ASSERT(NONULL(ele), "NULL as argument of function element_type_specific()\n");
    switch (ele->type) {
        case LINE:
            ele->dim = 1;
            ele->n_sides = 2;
            ele->n_nodes = 2;
            break;
        case TRIANGLE:
            ele->dim = 2;
            ele->n_sides = 3;
            ele->n_nodes = 3;
            break;
        case TETRAHEDRON:
            ele->dim = 3;
            ele->n_sides = 4;
            ele->n_nodes = 4;
            break;
        default:
            xprintf(UsrErr, "Element %d is of the unsupported type %d\n", ele.id(), ele->type);
    }
}

/**
 * ALLOCATION OF ARRAYS IN STRUCT ELEMENT WHOSE SIZE DEPENDS ONLY ON ELEMENT'S
 * TYPE, NOT ON TOPOLOGY OF THE MESH
 */
void GmshMeshReader::element_allocation_independent(ElementFullIter ele) {
    int si, ni;

    ASSERT(NONULL(ele), "NULL as argument of function element_allocation_independent()\n");
    ele->node = (Node**) xmalloc(ele->n_nodes * sizeof (Node*));
    ele->side = (struct Side**) xmalloc(ele->n_sides * sizeof ( struct Side*));

    ele->rhs = (double*) xmalloc(ele->n_sides * sizeof ( double));
    ele->bas_alfa = (double *) xmalloc(ele->n_sides * sizeof ( double));
    ele->bas_beta = (double *) xmalloc(ele->n_sides * sizeof ( double));
    ele->bas_gama = (double *) xmalloc(ele->n_sides * sizeof ( double));
    ele->bas_delta = (double *) xmalloc(ele->n_sides * sizeof ( double));

    FOR_ELEMENT_NODES(ele, ni) {
        ele->node[ ni ] = NULL;
    }

    FOR_ELEMENT_SIDES(ele, si) {
        ele->side[ si ] = NULL;
        ele->rhs[ si ] = 0.0;
        ele->bas_alfa[ si ] = 0.0;
        ele->bas_beta[ si ] = 0.0;
        ele->bas_gama[ si ] = 0.0;
        ele->bas_delta[ si ] = 0.0;
    }
}

/**
 * private method for reading of elements - in process of implementation
 */
void GmshMeshReader::read_elements(FILE* in, Mesh * mesh) {
    xprintf(Msg, " - Reading elements...");

    char line[ LINE_SIZE ];

    skip_to(in, "$Elements");
    xfgets(line, LINE_SIZE - 2, in);
    int numElements = atoi(xstrtok(line));
    ASSERT(!(numElements < 1), "Number of elements < 1 in function read_node_list()\n");

    mesh->element.reserve(numElements);
    //init_element_list( mesh );

    for (int i = 0; i < numElements; ++i) {
        xfgets(line, LINE_SIZE - 2, in);
        parse_element_line(mesh->element, line, mesh);
    }
    ASSERT((mesh->n_elements() == numElements),
            "Number of created elements %d does not match number of elements %d in the input file.\n",
            mesh->n_elements(),
            numElements);

    xprintf(Msg, " %d elements readed. ", mesh->n_elements())/*orig verb 4*/;
    xprintf(Msg, "O.K.\n")/*orig verb 2*/;
}
