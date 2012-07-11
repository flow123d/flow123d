/*
 * f123dPrcsng.cc
 *
 *  Created on: Nov 19, 2010
 *      Author: Honza + Dalibor
 */

#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <fstream>
#include <iomanip>
#include <cmath>

#include <malloc.h>
#include <string.h>
#include <stdlib.h>

#include "yaml/yamlReader.hh"

#include "mesh/Physical.hh"
#include "mesh/Element.hh"
#include "mesh/EPoint.hh"
#include "mesh/ELine.hh"
#include "mesh/ETriangle.hh"
#include "mesh/ETetrahedron.hh"

#include "boundary/AbstractBoundary.hh"
#include "boundary/BoundaryData.hh"
#include "boundary/Boundary.hh"

using namespace std;

#define MAX_NODES_IN_ELEMENT 4
#define MAX_TAGS_IN_ELEMENT 3

void printFooter() {
    printf("********************************************************************************\n");
    printf("***                           f123dTPeg -      End                           ***\n");
    printf("********************************************************************************\n");
}

void printHeader(int argc, char** argv) {
    /*
     * kontrola souboru s nastavenim
     */
    if (argc != 2) {
        printf("usage: %s <initFile.yml>\n", argv[0]);
        exit(1);
    }

    printf("********************************************************************************\n");
    printf("***                           f123dTPeg -    Start                           ***\n");
    printf("********************************************************************************\n");
}

void programExit(int errorLevel, int lineNum) {
    printf("\n");
    printf("--------------------------------------------------------------------------------\n");
    printf("   Error level: %d\n", errorLevel);
    printf("   Error line:  %d\n", lineNum);
    printFooter();

    exit(errorLevel);
}

void programExit() {
    printFooter();
    exit(0);
}

void readPhysical(ifstream* ifMesh, list<Physical*>* listPhysical, bool* changePhysicalId) {
    char fileLine[256];
    // -------------------------------------------------------------------------
    ifMesh->getline(fileLine, 256);
    if (strcmp(fileLine, "$PhysicalNames") == 0) {
        list<Physical*>* tmpListPhysical = new list<Physical*>();
        // ---------------------------------------------------------------------
        int numPhysical;
        *ifMesh >> numPhysical;
        ifMesh->getline(fileLine, 256);
        // ---------------------------------------------------------------------
        for (int i = 0; i < numPhysical; ++i) {
            int dim;
            *ifMesh >> dim;
            // -----------------------------------------------------------------
            int id;
            *ifMesh >> id;
            // -----------------------------------------------------------------
            ifMesh->getline(fileLine, 256);
            // -----------------------------------------------------------------
            // kontrola unikatnosti cislovani Physical
            for (list<Physical*>::const_iterator tmpPhysical = tmpListPhysical->begin(); tmpPhysical != tmpListPhysical->end(); tmpPhysical++) {
                if ((*tmpPhysical)->getId() == id) {
                    *changePhysicalId = true;
                    break;
                }
            }
            // -----------------------------------------------------------------
            tmpListPhysical->insert(tmpListPhysical->end(), new Physical(dim, id, fileLine));
        }
        // ---------------------------------------------------------------------
        // Kopirovani vrstev z pracovni mapy do listu.
        // V pripade, ze vrstvy nejsou unikatne cislovany, zajisti precislovani.
        for (list<Physical*>::const_iterator tmpPhysical = tmpListPhysical->begin(); tmpPhysical != tmpListPhysical->end(); tmpPhysical++) {
            int dim = (*tmpPhysical)->getDim();
            int id = (*tmpPhysical)->getId();
            if (*changePhysicalId) {
                id += 1000 * dim;
            }
            char* name = (*tmpPhysical)->getName();
            listPhysical->insert(listPhysical->end(), new Physical(dim, id, name));
        }
        // ---------------------------------------------------------------------
        ifMesh->getline(fileLine, 256);
        if (strcmp(fileLine, "$EndPhysicalNames") != 0) {
            cout << "Error reading mesh file      - found '" << fileLine << "' and wait for '$EndPhysicalNames'" << endl;
            programExit(1, __LINE__);
        }
    }
}

void readNodes(ifstream* ifMesh, map<int, Node*>* mapNode) {
    char fileLine[256];
    // -------------------------------------------------------------------------
    ifMesh->getline(fileLine, 256);
    if (strcmp(fileLine, "$Nodes") != 0) {
        cout << "Error reading mesh file " << endl;
        programExit(1, __LINE__);
    }
    // -------------------------------------------------------------------------
    int numNodes;
    *ifMesh >> numNodes;
    ifMesh->getline(fileLine, 256); // odstraneni zbytku radku
    cout << "   Pocet nodu: " << numNodes << endl;
    // -------------------------------------------------------------------------
    int nLabel;
    double x, y, z;
    for (int i = 0; i < numNodes; i++) { // nacitani vsech nodu
        *ifMesh >> nLabel;
        // ---------------------------------------------------------------------
        *ifMesh >> x;
        *ifMesh >> y;
        *ifMesh >> z;
        // ---------------------------------------------------------------------
        Node* node = new Node(nLabel, x, y, z);
        // ---------------------------------------------------------------------
        mapNode->insert(pair<int, Node*>(nLabel, node));
    }
    ifMesh->getline(fileLine, 256); // docteni zbytku radku
    // -------------------------------------------------------------------------
    ifMesh->getline(fileLine, 256);
    if (strcmp(fileLine, "$EndNodes") != 0) {
        cout << "Error reading mesh file      - found '" << fileLine << "' and wait for '$EndNodes'" << endl;
        programExit(1, __LINE__);
    }
}

void readElements(ifstream* ifMesh, map<int, Element*>* mapElement, map<int, Node*>* mapNode, bool* changePhysicalId) {
    char fileLine[256];
    // -------------------------------------------------------------------------
    int numPoints = 0;
    int numLines = 0;
    int numTriangles = 0;
    int numTetrahedrons = 0;
    // -------------------------------------------------------------------------
    ifMesh->getline(fileLine, 256);
    if (strcmp(fileLine, "$Elements") != 0) {
        cout << "Error reading mesh file " << endl;
        programExit(1, __LINE__);
    }
    // -------------------------------------------------------------------------
    int numElements;
    *ifMesh >> numElements;
    ifMesh->getline(fileLine, 256); // odstraneni zbytku radku
    cout << "   Pocet elementu: " << numElements << endl;
    // -------------------------------------------------------------------------
    int* tags = (int*) malloc(MAX_TAGS_IN_ELEMENT * sizeof (int));
    Node** nodes = (Node**) malloc(MAX_NODES_IN_ELEMENT * sizeof (Node*));
    // =========================================================================
    for (int ie = 0; ie < numElements; ie++) { // nacitani vsech elementu
        int eLabel;
        *ifMesh >> eLabel;
        // ---------------------------------------------------------------------
        short elmType;
        *ifMesh >> elmType;
        // ---------------------------------------------------------------------
        short numTags;
        *ifMesh >> numTags;
        // ---------------------------------------------------------------------
        for (int k = 0; k < numTags; k++) {
            int tag;
            *ifMesh >> tag;
            // -----------------------------------------------------------------
            *(tags + k) = tag;
        }
        // ---------------------------------------------------------------------
        for (int iin = 0; iin < Element::getNumNodes(elmType); iin++) {
            int nodeLabel;
            *ifMesh >> nodeLabel;
            // -----------------------------------------------------------------
            *(nodes + iin) = mapNode->find(nodeLabel)->second;
        }
        // =====================================================================
        // vytvoreni Elementu
        Element* element = NULL;
        // ---------------------------------------------------------------------
        switch (elmType) {
            case ELEMENT_TYPE_LINE:
                ++numLines;
                if (*changePhysicalId) {
                    *(tags) += 1000;
                }
                element = new ELine(eLabel, numTags, tags, nodes);
                break;
            case ELEMENT_TYPE_TRIANGLE:
                ++numTriangles;
                if (*changePhysicalId) {
                    *(tags) += 2000;
                }
                element = new ETriangle(eLabel, numTags, tags, nodes);
                break;
            case ELEMENT_TYPE_TETRAHEDRON:
                ++numTetrahedrons;
                if (*changePhysicalId) {
                    *(tags) += 3000;
                }
                element = new ETetrahedron(eLabel, numTags, tags, nodes);
                break;
            case ELEMENT_TYPE_POINT:
                ++numPoints;
                element = new EPoint(eLabel, numTags, tags, nodes);
                break;
        }
        // ---------------------------------------------------------------------
        mapElement->insert(pair<int, Element*>(eLabel, element));
    }
    // -------------------------------------------------------------------------
    free(nodes);
    free(tags);
    // -------------------------------------------------------------------------
    ifMesh->getline(fileLine, 256); // odstraneni zbytku radku
    // -------------------------------------------------------------------------
    ifMesh->getline(fileLine, 256);
    if (strcmp(fileLine, "$EndElements") != 0) {
        cout << "Error reading mesh file " << endl;
        programExit(1, __LINE__);
    }
    // -------------------------------------------------------------------------
    cout << "      num Points : " << numPoints << endl;
    cout << "      num Lines : " << numLines << endl;
    cout << "      num Triangles : " << numTriangles << endl;
    cout << "      num Tetrahedrons : " << numTetrahedrons << endl;
}

// *****************************************************************************
// ***   readMesh                                                            ***
// *****************************************************************************

void readMesh(const string fileName,
        char* meshFormat,
        list<Physical*>* listPhysical,
        bool* changePhysicalId,
        map<int, Node*>* mapNode,
        map<int, Element*>* mapElement) {
    // -------------------------------------------------------------------------
    cout << "Read mesh - Begin" << endl;
    // -------------------------------------------------------------------------
    char fileLine[256];
    // -------------------------------------------------------------------------
    // open input stream
    ifstream ifMesh(fileName.c_str(), ifstream::in); // jmeno souboru se siti
    if (!ifMesh.is_open()) {
        cout << "Error opening file " << fileName << endl;
        programExit(1, __LINE__);
    }
    // =========================================================================
    // read MESH FORMAT
    ifMesh.getline(fileLine, 256);
    if (strcmp(fileLine, "$MeshFormat") != 0) {
        cout << "Error reading mesh file " << endl;
        programExit(1, __LINE__);
    }
    // -------------------------------------------------------------------------
    ifMesh.getline(fileLine, 256);
    cout << "   GMSH format: " << fileLine << endl;
    strcpy(meshFormat, fileLine);
    //meshFormat[strlen(fileLine)] = '\0';
    // -------------------------------------------------------------------------
    ifMesh.getline(fileLine, 256);
    if (strcmp(fileLine, "$EndMeshFormat") != 0) {
        cout << "Error reading mesh file " << endl;
        programExit(1, __LINE__);
    }
    // =========================================================================
    readPhysical(&ifMesh, listPhysical, changePhysicalId);
    readNodes(&ifMesh, mapNode);
    readElements(&ifMesh, mapElement, mapNode, changePhysicalId);
    // -------------------------------------------------------------------------
    ifMesh.close();
    // -------------------------------------------------------------------------
    cout << "Read mesh - End" << endl;
}

double groupSize(map<int, list<int>* >* mapGroup, int groupLabel, map<int, Element*>* mapElement) {
    double sumSize = 0.0;
    list<int>* elmList = mapGroup->find(groupLabel)->second;
    for (list<int>::const_iterator iterElmLabel = elmList->begin(); iterElmLabel != elmList->end(); iterElmLabel++) {
        int elmLabel = (*iterElmLabel);
        Element* element = mapElement->find(elmLabel)->second;
        double size = element->getSize();
        sumSize += size;
    }
    return sumSize;
}

void writeTic(const string fileName, map<int, Element*>* mapElement, const vector<int>* excludePhysicalIds, const int numExcludedElm) {
    // ---------------------------------------------------------------------
    cout << "Write TIC: ";
    // ---------------------------------------------------------------------
    ofstream ofTic(fileName.c_str(), ofstream::out);
    if (!ofTic.is_open()) {
        cout << "Error opening file " << fileName << endl;
        programExit(1, __LINE__);
    }
    // ---------------------------------------------------------------------
    // zapis pocatecni podminky - hlavicka souboru
    ofTic << "$ConcentrationFormat" << endl << "1.0     0       8" << endl << "$EndConcentrationFormat" << endl << "$Concentrations" << endl;
    // ---------------------------------------------------------------------
    // zapis poctu elementu
    ofTic << mapElement->size() - numExcludedElm << endl;
    // ---------------------------------------------------------------------
    int ticCounter = 0;
    // ---------------------------------------------------------------------
    for (map<int, Element*>::const_iterator itElement = mapElement->begin(); itElement != mapElement->end(); itElement++) {
        Element* element = itElement->second;

        bool write = true;
        for (int j = 0; j < excludePhysicalIds->size(); j++) { // zjisteni jestli se ma element zapisovat do vysledku .msh
            if (element->getTag(0) == excludePhysicalIds->at(j)) {
                write = false;
                break;
            }
        }
        // -----------------------------------------------------------------
        if (write) {
            ofTic << ++ticCounter << "\t" << element->getLabel() << "\t" << "0.0" << endl;
        }
    }
    // ---------------------------------------------------------------------
    ofTic << "$EndConcentrations";
    // ---------------------------------------------------------------------
    ofTic.close();
    // ---------------------------------------------------------------------
    cout << "OK" << endl;
}

void readMaterial(const string fileName, map<int, map<char*, double>* >* mapMaterial) {
    const char* delim = " \t";
    char fileLine[256];
    bool sw = false;
    // -------------------------------------------------------------------------
    cout << "Read MTR - Begin" << endl;
    // ---------------------------------------------------------------------
    ifstream ifMatr(fileName.c_str(), ifstream::in); // jmeno souboru s materialovymi parametry
    if (!ifMatr.is_open()) {
        cout << "   Error opening file " << fileName << endl;
        programExit(1, __LINE__);
    }
    // ---------------------------------------------------------------------
    for (;;) {
        ifMatr.getline(fileLine, 256);
        // ---------------------------------------------------------------------
        if (strlen(fileLine) == 0) {
            break;
        }
        // ---------------------------------------------------------------------
        if (strcmp(fileLine, "$Geometry") == 0) {
            sw = true;
            continue;
        }
        // -----------------------------------------------------------------
        if (strcmp(fileLine, "$EndGeometry") == 0) {
            break;
        }
        // -----------------------------------------------------------------
        if (sw) {
            char* tokGroupLabel = strtok(fileLine, delim);
            char* tokGeoType = strtok(NULL, delim);
            char* tokValue = strtok(NULL, delim);
            // -------------------------------------------------------------
            int groupLabel = atoi(tokGroupLabel);
            double value = atof(tokValue);
            // -------------------------------------------------------------
            map<char*, double>* paramMap = mapMaterial->find(groupLabel)->second;
            if (paramMap == NULL) {
                paramMap = new map<char*, double>();
                mapMaterial->insert(pair<int, map<char*, double>* >(groupLabel, paramMap));
            } else {
                map<char*, double>::const_iterator iterParam = paramMap->find("Geometry");
                if (iterParam != paramMap->end()) {
                    cout << "   Multiply definition of Geometry parameter for group " << groupLabel << endl;
                    programExit(1, __LINE__);
                }
            }
            paramMap->insert(pair<char*, double>("Geometry", value));
        }
    }
    // ---------------------------------------------------------------------
    ifMatr.close();
    // ---------------------------------------------------------------------
    cout << "   Read material parameters for " << mapMaterial->size() << " groups" << endl;
    // ---------------------------------------------------------------------
    cout << "Read MTR - End" << endl;
}

//void writeMesh(const string fileName, const vector<int>* excludePhysicalIds, const int numExcludedElm) {

void writeMesh(const string fileName,
        char* meshFormat,
        list<Physical*>* listPhysical,
        map<int, Node*>* mapNode,
        map<int, Element*>* mapElement,
        const vector<int>* excludePhysicalIds,
        const int numExcludedElm) {
    cout << "Write MESH: ";
    // open output stream
    ofstream ofMesh(fileName.c_str(), ofstream::out);
    if (!ofMesh.is_open()) {
        cout << "Error opening file " << fileName << endl;
        programExit(1, __LINE__);
    }
    // =========================================================================
    // write MESH FORMAT
    if (meshFormat != NULL) {
        ofMesh << "$MeshFormat" << endl;
        ofMesh << meshFormat << endl;
        ofMesh << "$EndMeshFormat" << endl;
    }
    // =========================================================================
    // write PHYSICAL
    if (listPhysical->size() != 0) {
        ofMesh << "$PhysicalNames" << endl;
        ofMesh << listPhysical->size() << endl;
        for (list<Physical*>::const_iterator physical = listPhysical->begin(); physical != listPhysical->end(); physical++) {
            // ofMesh << (*iterPhysical)->dim << " " << ((*iterPhysical)->id * 100) << " " << (*iterPhysical)->fileName << endl;
            ofMesh << (*physical)->getDim() << " " << (*physical)->getId() << " " << (*physical)->getName() << endl;
        }
        ofMesh << "$EndPhysicalNames" << endl;
    }
    // =========================================================================
    // write NODES
    ofMesh << "$Nodes" << endl;
    ofMesh << mapNode->size() << endl;
    // -------------------------------------------------------------------------
    for (map<int, Node*>::const_iterator itNode = mapNode->begin(); itNode != mapNode->end(); itNode++) { // prochazeni vsech nodu
        Node* node = itNode->second;
        ofMesh << setprecision(15) << node->getLabel() << " " << node->getCoor(COOR_X) << " " << node->getCoor(COOR_Y) << " " << node->getCoor(COOR_Z) << endl;
    }
    // -------------------------------------------------------------------------
    ofMesh << "$EndNodes" << endl;
    // =========================================================================
    ofMesh << "$Elements" << endl;
    // -------------------------------------------------------------------------
    // zapis poctu elementu
    ofMesh << mapElement->size() - numExcludedElm << endl;
    // -------------------------------------------------------------------------
    for (map<int, Element*>::const_iterator itElement = mapElement->begin(); itElement != mapElement->end(); itElement++) {
        Element* element = itElement->second;

        bool write = true;
        for (int j = 0; j < excludePhysicalIds->size(); j++) { // zjisteni jestli se ma element zapisovat do vysledku .msh
            if (element->getTag(0) == excludePhysicalIds->at(j)) {
                write = false;
                break;
            }
        }
        // ---------------------------------------------------------------------
        if (write) {
            ofMesh << element->getLabel() << " ";
            ofMesh << element->getElementType() << " ";
            ofMesh << element->getNumTags() << " ";
            for (int iiTag = 0; iiTag < element->getNumTags(); iiTag++) {
                ofMesh << element->getTag(iiTag) << " ";
            }
            for (int iin = 0; iin < element->getNumNodes(); ++iin) {
                ofMesh << element->getNode(iin)->getLabel() << " ";
            }
            // -----------------------------------------------------------------
            ofMesh << endl;
        }
    }
    // -------------------------------------------------------------------------
    ofMesh << "$EndElements" << endl;
    // -------------------------------------------------------------------------
    ofMesh.close();
    // -------------------------------------------------------------------------
    cout << "OK" << endl;
}

void readNgh(const string fileName, map<int, list< int* >* >* mapElmNgh) {
    cout << "Read NGH - Begin" << endl;
    // -------------------------------------------------------------------------
    char fileLine[256];
    // -------------------------------------------------------------------------
    ifstream ifNgh(fileName.c_str(), ifstream::in); // jmeno souboru se siti
    if (!ifNgh.is_open()) {
        cout << "   Error opening file " << fileName << endl;
        programExit(1, __LINE__);
    }
    // -------------------------------------------------------------------------
    for (int i = 0; i < 4; i++) {
        ifNgh.getline(fileLine, 256); // odstraneni prvnich 5 radku souboru sousednosti
    }
    // -------------------------------------------------------------------------
    int numNeighbours;
    ifNgh >> numNeighbours;
    // -------------------------------------------------------------------------
    int numNGH20 = 0;
    for (int i = 0; i < numNeighbours; i++) { // prochazeni celeho souboru sosednosti
        // bude nas zajimat pouze sousednost 20, ostatni zahazujeme
        int nghId;
        ifNgh >> nghId;
        // ---------------------------------------------------------------------
        short type;
        ifNgh >> type;
        // ---------------------------------------------------------------------
        if (type != 20) {
            char fileLine[256];
            ifNgh.getline(fileLine, 256); // prectu zbytek radku a vyhodim jej
            continue;
        }
        // ---------------------------------------------------------------------
        ++numNGH20;
        // ---------------------------------------------------------------------
        int eid1;
        ifNgh >> eid1;
        // ---------------------------------------------------------------------
        int eid2;
        ifNgh >> eid2;
        // ---------------------------------------------------------------------
        int sid2;
        ifNgh >> sid2;
        // ---------------------------------------------------------------------
        ifNgh.getline(fileLine, 256); // prectu zbytek radku a vyhodim jej
        // =====================================================================
        list< int* >* listList = mapElmNgh->find(eid1)->second;
        if (listList == NULL) {
            int* eid_sid = (int*) malloc(2 * sizeof (int));
            *(eid_sid + 0) = eid2;
            *(eid_sid + 1) = sid2;
            listList = new list< int* >(1, eid_sid);
            mapElmNgh->insert(pair<int, list< int* >*>(eid1, listList));
        } else {
            int* eid_sid = (int*) malloc(2 * sizeof (int));
            *(eid_sid + 0) = eid2;
            *(eid_sid + 1) = sid2;
            listList->insert(listList->end(), eid_sid);
        }
    }
    // -------------------------------------------------------------------------
    ifNgh.close();
    // -------------------------------------------------------------------------
    cout << "   Read neighbours for " << mapElmNgh->size() << " elements" << endl;
    // ---------------------------------------------------------------------
    cout << "Read NGH - End" << endl;
}

void writeFBc(const string fileName, list<Boundary*>* listBoundary) {
    cout << "Write FBC: ";
    // -------------------------------------------------------------------------
    ofstream ofBcd(fileName.c_str(), ofstream::out); // jmeno souboru se siti
    if (!ofBcd.is_open()) {
        cout << "Error opening output file " << fileName.c_str() << endl;
        programExit(1, __LINE__);
    }
    // -------------------------------------------------------------------------
    ofBcd << "$BoundaryFormat" << endl;
    ofBcd << "1.0  0  8" << endl;
    ofBcd << "$EndBoundaryFormat" << endl;
    ofBcd << "$BoundaryConditions" << endl;
    // -------------------------------------------------------------------------
    // ofBcd << dirichletBC.size() + bcdNeumannCounter2 << endl;
    ofBcd << listBoundary->size() << endl;
    // -------------------------------------------------------------------------
    // vypsat okrajove podminky do souboru
    for (list<Boundary*>::const_iterator iterBoundary = listBoundary->begin(); iterBoundary != listBoundary->end(); iterBoundary++) {
        int bcType = (*iterBoundary)->getType();
        // ---------------------------------------------------------------------
        ofBcd << (*iterBoundary)->getId() << " ";
        ofBcd << (*iterBoundary)->getType() << " ";
        ofBcd << (*iterBoundary)->getValue() << " ";
        if ((*iterBoundary)->getType() == NEWTON_BC) {
            ofBcd << (*iterBoundary)->getSigma() << " ";
        }
        ofBcd << (*iterBoundary)->getWhere() << " ";
        ofBcd << (*iterBoundary)->getElmId() << " ";
        ofBcd << (*iterBoundary)->getSidId() << " ";
        ofBcd << "1" << " "; // number of Tags
        ofBcd << (*iterBoundary)->getTag() << endl;
    }
    // -------------------------------------------------------------------------
    ofBcd << "$EndBoundaryConditions" << endl;
    ofBcd.close();
    // -------------------------------------------------------------------------
    cout << "OK" << endl;
}

void writeTBc(const string fileName) {
    cout << "Write TBC: ";
    // -------------------------------------------------------------------------
    ofstream ofTbc(fileName.c_str(), ofstream::out); // jmeno souboru se siti
    if (!ofTbc.is_open()) {
        cout << "Error opening output file " << fileName.c_str() << endl;
        programExit(1, __LINE__);
    }
    // -------------------------------------------------------------------------
    ofTbc << "$Transport_BCDFormat" << endl;
    ofTbc << "1.0  0  8" << endl;
    ofTbc << "$EndTransport_BCDFormat" << endl;
    ofTbc << "$Transport_BCD" << endl;
    // -------------------------------------------------------------------------
    // ofTbc << listBoundary->size() + bcdNeumannCounter2 << endl;
    // -------------------------------------------------------------------------
    // vypsat okrajove podminky do souboru
    int transportBcdCounter = 0;
    // for (vector<Dirichlet>::const_iterator iterDirichlet = dirichletBC.begin(); iterDirichlet < dirichletBC.end(); iterDirichlet++) { // zapiseme vsechny okrajove podminky dirichletova typu
    // ofTbc << transportBcdCounter++ << "\t" << (*iterDirichlet).id << "\t";
    // ofTbc.setf(ios::fixed, ios::floatfield);
    // ofTbc.precision(5);
    // ofTbc << 0.0 << endl;
    // }
    // -------------------------------------------------------------------------
    double f = 1.00;
    // for (vector<Neumann>::const_iterator iterNeumann = neumannV2.begin(); iterNeumann < neumannV2.end(); iterNeumann++) {
    // ofTbc << transportBcdCounter++ << "\t" << (*iterNeumann).id + dirichletBC.size() << "\t";
    // ofTbc.setf(ios::fixed, ios::floatfield);
    // ofTbc.precision(5);
    // ofTbc << f << endl;
    // }
    // -------------------------------------------------------------------------
    ofTbc << "$EndTransport_BCD" << endl;
    ofTbc.close();
    // -------------------------------------------------------------------------
    cout << "OK" << endl;
}

void readTerrain(const string fileName, list< vector<double>* >* terenList) {
    cout << "Read terrain: ";
    // -------------------------------------------------------------------------
    char fileLine[256];
    // -------------------------------------------------------------------------
    ifstream ifTeren(fileName.c_str(), ifstream::in); // jmeno souboru s reliefem
    if (!ifTeren.is_open()) {
        cout << "Error opening file '" << fileName << "'" << endl;
        programExit(1, __LINE__);
    }
    // -------------------------------------------------------------------------
    ifTeren.getline(fileLine, 256); // odstraneni prvniho radku s popiskama
    // -------------------------------------------------------------------------
    char semicolon;
    // -------------------------------------------------------------------------
    while (ifTeren.good()) {
        ifTeren >> semicolon;
        if (semicolon == '\n') {
            continue;
        }
        // ---------------------------------------------------------------------
        int id;
        ifTeren >> id;
        // ---------------------------------------------------------------------
        ifTeren >> semicolon;
        // ---------------------------------------------------------------------
        double x;
        ifTeren >> x;
        // ---------------------------------------------------------------------
        ifTeren >> semicolon;
        // ---------------------------------------------------------------------
        double y;
        ifTeren >> y;
        // ---------------------------------------------------------------------
        ifTeren >> semicolon;
        // ---------------------------------------------------------------------
        double z;
        ifTeren >> z;
        // ---------------------------------------------------------------------
        vector<double>* vectorTmp = new vector<double>(3);
        vectorTmp->assign(0, x);
        vectorTmp->assign(1, y);
        vectorTmp->assign(2, z);
        // ---------------------------------------------------------------------
        terenList->insert(terenList->end(), vectorTmp);
        // ---------------------------------------------------------------------
        ifTeren.getline(fileLine, 256); // odstraneni zbytku radku
    }
    // -------------------------------------------------------------------------
    ifTeren.close();
    // -------------------------------------------------------------------------
    cout << "OK     Pocet bodu terenu: " << terenList->size() << endl;
}

void readResults(const string fileName) {
    cout << "Read results - Begin" << endl;
    cout << "   NOT completly implemented" << endl;
    return;
    // -------------------------------------------------------------------------
    char fileLine[256];
    // -------------------------------------------------------------------------
    ifstream ifResults(fileName.c_str(), ifstream::in); // jmeno souboru se siti
    if (!ifResults.is_open()) {
        cout << "FlowResults file '" << fileName << "' not found" << endl;
        programExit(1, __LINE__);
    }
    // -------------------------------------------------------------------------
    for (int i = 0; i < 5; i++) {
        ifResults.getline(fileLine, 256); // odstraneni prvnich 5 radku souboru vysledku
    }
    // -------------------------------------------------------------------------
    int numResults;
    ifResults >> numResults;
    ifResults.getline(fileLine, 256); // odstraneni zbytku radku
    // -------------------------------------------------------------------------
    for (int i = 0; i < numResults; i++) {
        int cit;
        ifResults >> cit; // necteni cisla vysledku
        // ---------------------------------------------------------------------
        int elmId;
        ifResults >> elmId; // necteni Id elementu
        // ---------------------------------------------------------------------
        short elmNsides;
        ifResults >> elmNsides; // necteni poctu sten daneho elementu
        // ---------------------------------------------------------------------
        double clear;
        ifResults >> clear; // neni potreba, vyhodit
        // ---------------------------------------------------------------------
        for (int j = 0; j < elmNsides; j++) ifResults >> clear; // odmazani skalaru na stenach
        vector<double> sidFlux(elmNsides, 0.0);
        for (int j = 0; j < elmNsides; j++) {
            ifResults >> sidFlux[j];
        }
        // ---------------------------------------------------------------------
        for (int j = 0; j < 3; j++) ifResults >> clear; // odmazani vectoru na elementu
        // ---------------------------------------------------------------------
        int clearI;
        ifResults >> clearI;
        // ---------------------------------------------------------------------
        // for (vector<Neumann>::const_iterator neumannFlux = neumannV2.begin(); neumannFlux < neumannV2.end(); neumannFlux++) {
        // if ((*neumannFlux).elmId == elmId) {
        // (*neumannFlux).flux = sidFlux[(*neumannFlux).sidId] / 1e4;
        // }
        // }
    }
    // -------------------------------------------------------------------------
    ifResults.close();
    // -------------------------------------------------------------------------
    cout << "Read results - End" << endl;
}

void debugWriteMesh(map<int, Element*>* mapElement) {
    cout << "DEBUG: Write mesh - End" << endl;
    // -------------------------------------------------------------------------
    int ie = 0;
    // -------------------------------------------------------------------------
    for (map<int, Element*>::const_iterator itElement = mapElement->begin(); itElement != mapElement->end(); itElement++) {
        Element* element = itElement->second;
        int numNodes = element->getNumNodes();
        // ---------------------------------------------------------------------
        if (((ie % 10000) == 0) || (ie < 10)) {
            cout << "   Label: " << element->getLabel();
            cout << "     tagP: " << element->getTag(PHYSICAL_GROUP);
            cout << "     tagE: " << element->getTag(ELEMENT_GROUP);
            cout << "     nn: " << numNodes;
            for (int iin = 0; iin < numNodes; ++iin) {
                cout << "     n" << (iin + 1) << ": " << element->getNode(iin)->getLabel();
            }
            cout << endl;
        }
        // ---------------------------------------------------------------------
        ++ie;
    }
    // -------------------------------------------------------------------------
    cout << "DEBUG: Write mesh - End" << endl;
}

void debugWriteNgh(map<int, list< int* >* >* mapElmNgh) {
    cout << "DEBUG: Write NGH - Begin" << endl;
    // -------------------------------------------------------------------------
    for (map<int, list< int* >* >::const_iterator iterENL = mapElmNgh->begin(); iterENL != mapElmNgh->end(); iterENL++) {
        list< int* >* listList = iterENL->second;
        // -------------------------------------------------------------------------
        if (listList->size() > 1) {
            cout << "   " << iterENL->first << " " << listList->size() << endl;
            // -------------------------------------------------------------------------
            for (list< int* >::const_iterator nList = listList->begin(); nList != listList->end(); nList++) {
                cout << "      " << *(*nList + 0) << "   " << *(*nList + 1) << endl;
            }
        }
    }
    // -------------------------------------------------------------------------
    cout << "DEBUG: Write Ngh - End" << endl;
}

void terrainNodeModification(map<int, Node*>* mapNode, list< vector<double>* >* listTerrain) {
    cout << "Node terrain modification - Begin";
    // -------------------------------------------------------------------------
    int modifNodes = 0;
    double eps = 45; // polomer kruznice ve ktere je nejblizsi bod ke vsem nodum v siti
    int bottomLayerZ = -600; // z souradnice spodni vrstvy
    int topLayerZ = 550; // z souradnice vrchni vrstvy
    // -------------------------------------------------------------------------
    for (map<int, Node*>::const_iterator itNode = mapNode->begin(); itNode != mapNode->end(); itNode++) { // prochazeni vsech nodu
        Node* node = itNode->second;
        list< vector<double>* >::const_iterator iterTeren = listTerrain->begin();
        double minNorm = 1e10;
        for (; iterTeren != listTerrain->end(); iterTeren++) {
            // cout << "x = " << setprecision(15) << (*iter)[0] << "; y = " << (*iter)[1] << "; z = " << (*iter)[2] << endl;
            double norm = 0.0;
            // minNorm = minNorm > (norm = abs(sqrt(pow(((*iterTeren)[0] - nodeList[i].x), 2) + pow(((*iterTeren)[1] - nodeList[i].y), 2)))) ? norm : minNorm;
            minNorm = minNorm > (norm = abs(sqrt(pow(((*iterTeren)->at(0) - node->getCoor(COOR_X)), 2) + pow(((*iterTeren)->at(1) - node->getCoor(COOR_Y)), 2)))) ? norm : minNorm;
            if (minNorm < eps) {
                ++modifNodes;
                break;
            }
        }
        node->setCoor(COOR_Z, ((node->getCoor(COOR_Z) - bottomLayerZ) / (((double) topLayerZ) - bottomLayerZ))*((*iterTeren)->at(2) - bottomLayerZ) + bottomLayerZ); // vypocteni nove z souradnice dle reliefu
        node->setLevel((*iterTeren)->at(2) - node->getCoor(COOR_Z));
    }
    // -------------------------------------------------------------------------
    cout << "    Modificate " << modifNodes << " nodes" << endl;
    // -------------------------------------------------------------------------
    cout << "Node terrain modification - End";
}

void countElementMassPoint(map<int, Element*>* mapElement) {
    cout << "Vypocet teziste elementu: ";
    // dopocet teziste elementu
    //    pro vsechy elementy
    for (map<int, Element*>::const_iterator itElement = mapElement->begin(); itElement != mapElement->end(); itElement++) {
        Element* element = itElement->second;

        vector<double> centroid(4, 0.0);
        for (int iin = 0; iin < element->getNumNodes(); iin++) {
            Node* node = element->getNode(iin);
            // ------------------------------------------------------------------
            centroid[0] += node->getCoor(COOR_X);
            centroid[1] += node->getCoor(COOR_Y);
            centroid[2] += node->getCoor(COOR_Z);
            centroid[3] += node->getLevel();
        }

        double devide = element->getNumNodes();

        for (int iCoor = 0; iCoor < 3; ++iCoor) {
            centroid[iCoor] /= devide;
        }
        element->T = new Node(-1, centroid[COOR_X], centroid[COOR_Y], centroid[COOR_Z]);

        element->T->setLevel(centroid[3] / devide);
    }
    cout << "OK" << endl;
}

void markFictiveElements(map<int, Element*>* mapElement, int excludeType, vector<int>* excludeIds, int* numExcluded) {
    cout << "Marking of fictive elements (type " << excludeType << "): ";
    // -------------------------------------------------------------------------
    int excludeCounter = 0;
    // -------------------------------------------------------------------------
    for (map<int, Element*>::const_iterator itElement = mapElement->begin(); itElement != mapElement->end(); itElement++) {
        Element* element = itElement->second;
        // ---------------------------------------------------------------------
        if (element->isFictive()) {
            continue;
        }
        // ---------------------------------------------------------------------
        int tag = element->getTag(excludeType);
        for (int j = 0; j < excludeIds->size(); j++) {
            if (tag == excludeIds->at(j)) {
                excludeCounter++;
                element->setFictive(true);
                break;
            }
        }
    }
    // -------------------------------------------------------------------------
    *numExcluded += excludeCounter;
    // -------------------------------------------------------------------------
    cout << " marked " << excludeCounter << " elements" << endl;
}

void createGroupElements(map<int, Element*>* mapElement, int iTag, map<int, list<int>* >* mapGroup) {
    cout << "Create Group of elements (type " << iTag << "): ";
    // -------------------------------------------------------------------------
    for (map<int, Element*>::const_iterator itElement = mapElement->begin(); itElement != mapElement->end(); itElement++) {
        Element* element = itElement->second;
        // ---------------------------------------------------------------------
        int elmLabel = element->getLabel();
        int groupLabel = element->getTag(iTag);
        // ---------------------------------------------------------------------
        map<int, list<int>*>::const_iterator iterListElm = mapGroup->find(groupLabel);
        list<int>* listElm = iterListElm->second;
        // ---------------------------------------------------------------------
        if (iterListElm == mapGroup->end()) {
            listElm = new list<int>();
            mapGroup->insert(pair<int, list<int>* >(groupLabel, listElm));
        }
        // ---------------------------------------------------------------------
        listElm->insert(listElm->end(), elmLabel);
    }
    // -------------------------------------------------------------------------
    cout << "OK" << endl;
}

void writeGroupElements(map<int, list<int>* >* mapGroup, map<int, Element*>* mapElement) {
    cout << "Write Group of elements - Begin" << endl;
    // -------------------------------------------------------------------------
    int groupSum = 0;
    // ---------------------------------------------------------------------
    for (map<int, list<int>* >::const_iterator iterGroup = mapGroup->begin(); iterGroup != mapGroup->end(); iterGroup++) {
        int groupLabel = iterGroup->first;
        list<int>* listElm = iterGroup->second;
        // -----------------------------------------------------------------
        cout << "   Group: " << groupLabel;
        cout << "   num elements: " << listElm->size();
        cout << "   group size: " << groupSize(mapGroup, groupLabel, mapElement);
        cout << endl;
        // -----------------------------------------------------------------
        groupSum += listElm->size();
        // -----------------------------------------------------------------
        if (false) {
            cout << "      DEBUG: ";
            if (listElm->size() < 20) {
                for (list<int>::const_iterator iterElm = listElm->begin(); iterElm != listElm->end(); iterElm++) {
                    int elmLabel = (*iterElm);
                    cout << " " << elmLabel;
                }
                cout << endl;
            }
        }
    }
    // -------------------------------------------------------------------------
    cout << "   Group sum: " << groupSum << endl;
    // -------------------------------------------------------------------------
    cout << "Write Group of elements - End" << endl;
}

/**************************************************************************
 ***   Boundary Condition calculation                                     *
 **************************************************************************/
void createBoundaryConditions(list<BoundaryData*>* listBoundaryData,
        map<int, list<int>* >* mapGroup,
        map<int, Element*>* mapElement,
        map<int, map<char*, double>* >* mapMaterial,
        map<int, list< int* >* >* mapElmNgh,
        list<Boundary*>* listBoundary) {
    // -------------------------------------------------------------------------
    cout << "Boundary Condition calculation - Begin" << endl;
    // -------------------------------------------------------------------------
    int bcCounter = 0;
    // -------------------------------------------------------------------------
    int dirichletPHCounter = 0;
    int dirichletHPCounter = 0;
    int neumannPSCounter = 0;
    int neumannSumCounter = 0;
    int newtonPHCounter = 0;
    int newtonHPCounter = 0;
    // -------------------------------------------------------------------------
    for (list<BoundaryData*>::const_iterator iterBoundaryData = listBoundaryData->begin(); iterBoundaryData != listBoundaryData->end(); iterBoundaryData++) {
        int bcTag = (*iterBoundaryData)->getTag();
        int bcDataType = (*iterBoundaryData)->getType();
        double bcValue = (*iterBoundaryData)->getValue();
        double bcSigma = 0.0;
        // ---------------------------------------------------------------------
        switch (bcDataType) {
            case DIRICHLET_PIEZOMETRIC_HEAD:
                cout << "   Group with Dirichlet BC (piezometric head) - group:" << bcTag;
                cout << "   value:" << bcValue << endl;
                break;
            case DIRICHLET_HEAD_PRESSURE:
                cout << "   Group with Dirichlet BC (head pressure)    - group:" << bcTag;
                cout << "   value:" << bcValue << endl;
                break;
            case NEUMANN_PER_SQUARE:
                cout << "   Group with Neumann BC (pre square)         - group:" << bcTag;
                cout << "   value:" << bcValue << endl;
                break;
            case NEUMANN_SUM:
                cout << "   Group with Neumann BC (sum)                - group:" << bcTag;
                cout << "   value:" << bcValue << endl;
                break;
            case NEWTON_PIEZOMETRIC_HEAD:
                cout << "   Group with Newton BC (piezometric head)    - group:" << bcTag;
                bcSigma = (*iterBoundaryData)->getSigma();
                cout << "   value:" << bcValue << "   sigma:" << bcSigma << endl;
                break;
            case NEWTON_HEAD_PRESSURE:
                cout << "   Group with Newton BC (head pressure)       - group:" << bcTag;
                bcSigma = (*iterBoundaryData)->getSigma();
                cout << "   value:" << bcValue << "   sigma:" << bcSigma << endl;
                break;
        }
        // ---------------------------------------------------------------------
        list<int>* listElm = mapGroup->find(bcTag)->second;
        // ---------------------------------------------------------------------
        for (list<int>::const_iterator iterElm = listElm->begin(); iterElm != listElm->end(); iterElm++) {
            int elmLabel = (*iterElm);
            Element* element = mapElement->find(elmLabel)->second;
            // -----------------------------------------------------------------
            list< int* >* listElmNgh = mapElmNgh->find(elmLabel)->second;
            // -----------------------------------------------------------------
            if (listElmNgh == NULL) {
                cout << "WARNING:   element: " << elmLabel << " have ZERO neigbours" << endl;
                continue;
            }
            // -----------------------------------------------------------------
            for (list< int* >::const_iterator iterList = listElmNgh->begin(); iterList != listElmNgh->end(); iterList++) {
                int eid2 = *(*iterList + 0);
                int sid2 = *(*iterList + 1);
                // -------------------------------------------------------------
                Element* elmNbr = mapElement->find(eid2)->second;
                // -------------------------------------------------------------
                if (elmNbr->isFictive()) {
                    continue;
                }
                // -------------------------------------------------------------
                double value;
                int bcType = -1;
                // -------------------------------------------------------------
                // Special variables for Neumann SUM flow BC
                double geoValue;
                int groupElm2;
                map<char*, double>* mapParam;
                // -------------------------------------------------------------
                switch (bcDataType) {
                    case DIRICHLET_PIEZOMETRIC_HEAD:
                        ++dirichletPHCounter;
                        value = bcValue - element->T->getCoor(COOR_Z); // prepocet PIZEOMETRIC HEAD -> HEAD PRESSURE
                        bcType = DIRICHLET_BC;
                        break;
                    case DIRICHLET_HEAD_PRESSURE:
                        ++dirichletHPCounter;
                        value = bcValue; // HEAD PRESSURE
                        bcType = DIRICHLET_BC;
                        break;
                        // -----------------------------------------------------
                    case NEUMANN_PER_SQUARE:
                        ++neumannPSCounter;
                        value = bcValue; //* element->getSize(); // vypocet VALUE_PER_SQUARE
                        bcType = NEUMANN_BC;
                        break;
                    case NEUMANN_SUM:
                        ++neumannSumCounter;
                        if (mapMaterial != NULL) {
                            groupElm2 = mapElement->find(eid2)->second->getTag(PHYSICAL_GROUP);
                            mapParam = mapMaterial->find(groupElm2)->second;
                            if (mapParam == NULL) {
                                cout << "   WARNING: Neumann SUM boundary condition" << endl;
                                cout << "            For element group(" << groupElm2 << ") is NOT defined Geometry parameter in file *.mtr" << endl;
                                cout << "            Geometry parameter is set = 1.0" << endl;
                                geoValue = 1.0;
                            } else {
                                geoValue = mapParam->find("Geometry")->second;
                            }
                        } else {
                            cout << "   WARNING: Neumann SUM boundary condition" << endl;
                            cout << "            File *.mtr is NOT defined" << endl;
                            cout << "            Geometry parameter is set = 1.0" << endl;
                            geoValue = 1.0;
                        }
                        value = bcValue * (1.0 / (geoValue * groupSize(mapGroup, bcTag, mapElement))); // vypocet podilu Value pro dany Elm
                        bcType = NEUMANN_BC;
                        break;
                        // -----------------------------------------------------
                    case NEWTON_PIEZOMETRIC_HEAD:
                        ++newtonPHCounter;
                        value = bcValue - element->T->getCoor(COOR_Z); // prepocet PIZEOMETRIC HEAD -> HEAD PRESSURE
                        bcType = NEWTON_BC;
                        break;
                    case NEWTON_HEAD_PRESSURE:
                        ++newtonHPCounter;
                        value = bcValue; // HEAD PRESSURE
                        bcType = NEWTON_BC;
                        break;
                        // -----------------------------------------------------
                    default:
                        cout << "   WARNING - Unknown type of Boundary Condition: " << bcDataType << endl;
                        break;
                }
                // -------------------------------------------------------------
                int where = CONDITION_ON_SIDE;
                int elmId = eid2;
                int sidId = sid2;
                int tag = bcTag;
                // -------------------------------------------------------------
                listBoundary->insert(listBoundary->end(), new Boundary(bcCounter, bcType, value, bcSigma, where, elmId, sidId, tag));
                // -------------------------------------------------------------
                ++bcCounter;
            }
        }
    }
    // -------------------------------------------------------------------------
    cout << "--------------------------------------------" << endl;
    cout << "   Pocet Dirichlet BC PiezoHead: " << dirichletPHCounter << endl;
    cout << "   Pocet Dirichlet BC HeadPressure: " << dirichletHPCounter << endl;
    cout << "   Pocet Neumann BC pre Square: " << neumannPSCounter << endl;
    cout << "   Pocet Neumann BC Sum flow: " << neumannSumCounter << endl;
    cout << "   Pocet Newton BC PiezoHead: " << newtonPHCounter << endl;
    cout << "   Pocet Newton BC HeadPressure: " << newtonHPCounter << endl;
    cout << "--------------------------------------------" << endl;
    cout << "   Pocet BC celkem: " << listBoundary->size() << endl;
    // -------------------------------------------------------------------------
    cout << "Boundary Condition calculation - End" << endl;
}

int main(int argc, char** argv) {
    printHeader(argc, argv);

    // =========================================================================
    // =====                                                               =====
    // =========================================================================
    cout << "Read YML file - Begin" << endl;
    // -------------------------------------------------------------------------
    list<BoundaryData*>* listBoundaryData = new list<BoundaryData*>();
    // -------------------------------------------------------------------------
    /*
     * Vytvoreni objektu na cteni nastaveni
     */
    YamlReader yamlReader(argv[1]);
    cout << "      YmlReader(" << argv[1] << ") - created" << endl;
    // =========================================================================
    /* jmeno vstupniho souboru s terenem  - id: input > terrain > path */
    string terrainName;
    {
        string tmpPath[] = {"input", "terrain", "path"};
        cout << "      - reading 'terrain' - start" << endl;
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), &terrainName)) {
            cout << "   input: terrain <fileName> is not defined" << endl;
        }
    }
    cout << "      - reading 'terrain' - end" << endl;
    // -------------------------------------------------------------------------
    // vektor rozhodnych hodnot
    /*
     * vektor rozhodnych hodnot pro prenastaveni Physical Entity id
     * - pro moznost nastaveni materialovych parametru dle hloubky pod povrchem
     *
     * id: settings > decisiveDepth
     */
    vector<int>* devLayer = new vector<int>();
    {
        string tmpPath[] = {"settings", "decisiveDepth"};
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), devLayer)) {
            cout << "   settings: decisiveDepth are not defined" << endl;
        }
    }
    cout << "      - decisiveDepth" << endl;
    // -------------------------------------------------------------------------
    /* jmeno vstupniho souboru se siti  - id: input > mesh > path */
    string inputMesh;
    {
        string tmpPath[] = {"input", "mesh", "path"};
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), &inputMesh)) {
            cout << "   input: mesh <fileName> is not defined" << endl;
            programExit(1, __LINE__);
        }
    }
    // -------------------------------------------------------------------------
    /* jmeno vstupniho souboru s materialovymi mparametry  - id: input > matr > path */
    string inputMaterial;
    {
        string tmpPath[] = {"input", "matr", "path"};
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), &inputMaterial)) {
            cout << "   input: matr <fileName> is not defined" << endl;
        }
    }
    // -------------------------------------------------------------------------
    /* jmeno vystupniho souboru se siti  - id: output > mesh > path */
    string outputMesh;
    {
        string tmpPath[] = {"output", "mesh", "path"};
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), &outputMesh)) {
            cout << "   output: mesh: <fileName> is not defined" << endl;
            programExit(1, __LINE__);
        }
    }
    // -------------------------------------------------------------------------
    /* id odebiranych physical entit - id: settings > mesh > removeEntities > physical */
    vector<int>* excludePhysicalIds = new vector<int>();
    {
        string tmpPath[] = {"settings", "mesh", "removeEntities", "physical"};
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), excludePhysicalIds)) {
            cout << "   settings: mesh: removeEntities: [physical]  are not defined" << endl;
        }
    }
    // -------------------------------------------------------------------------
    /* id odebiranych elementary entit - id: settings > mesh > removeEntities > elementary */
    vector<int>* excludeElementaryIds = new vector<int>();
    {
        string tmpPath[] = {"settings", "mesh", "removeEntities", "elementary"};
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), excludeElementaryIds)) {
            cout << "   settings: mesh: removeEntities: [elementary]  are not defined" << endl;
        }
    }
    // -------------------------------------------------------------------------
    /* jmeno vystupniho souboru s pocatecni podminkou  - id: output > tic > path */
    string outputTic;
    {
        string tmpPath[] = {"output", "tic", "path"};
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), &outputTic)) {
            cout << "   output: tic <fileName> is not defined" << endl;
        }
    }
    // -------------------------------------------------------------------------
    string inputResults;
    {
        string tmpPath[] = {"input", "flowResults", "path"};
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), &inputResults)) {
            cout << "   input: flowResults <fileName> is not defined" << endl;
        }
    }
    // -------------------------------------------------------------------------
    /* jmeno vstupniho souboru se sousednosti  - id: input > ngh > path */
    string inputNgh;
    {
        string tmpPath[] = {"input", "ngh", "path"};
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), &inputNgh)) {
            cout << "   input: ngh <fileName> is not defined" << endl;
            programExit(1, __LINE__);
        }
    }
    // =========================================================================
    /* Dirichlet OKP - piezometric head - id: settings > dirichlet PiezometricHead */
    {
        string tmpPath[] = {"settings", "dirichlet", "piezometricHead"};
        vector<double>* dirichletPiezometricHead = new vector<double>();
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), dirichletPiezometricHead)) {
            cout << "   settings: dirichlet: piezometricHead [] are not defined" << endl;
        }

        for (vector<double>::const_iterator iterDirichlet = dirichletPiezometricHead->begin(); iterDirichlet < dirichletPiezometricHead->end(); iterDirichlet += 2) {
            int tag = (int) (*iterDirichlet);
            double value = (double) (*(iterDirichlet + 1));

            list<BoundaryData*>::iterator iter = listBoundaryData->end();
            int id = listBoundaryData->size();
            int type = DIRICHLET_PIEZOMETRIC_HEAD;

            listBoundaryData->insert(iter, new BoundaryData(id, type, value, tag));
        }
    }
    // -------------------------------------------------------------------------
    /* Dirichlet OKP - headPressure - id: settings > dirichlet HeadPressure */
    {
        string tmpPath[] = {"settings", "dirichlet", "headPressure"};
        vector<double>* dirichletHeadPressure = new vector<double>();
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), dirichletHeadPressure)) {
            cout << "   settings: dirichlet: headPressure [] are not defined" << endl;
        }

        for (vector<double>::const_iterator iterDirichlet = dirichletHeadPressure->begin(); iterDirichlet < dirichletHeadPressure->end(); iterDirichlet += 2) {
            int tag = (int) (*iterDirichlet);
            double value = (double) (*(iterDirichlet + 1));

            list<BoundaryData*>::iterator iter = listBoundaryData->end();
            int id = listBoundaryData->size();
            int type = DIRICHLET_HEAD_PRESSURE;

            listBoundaryData->insert(iter, new BoundaryData(id, type, value, tag));
        }
    }
    // -------------------------------------------------------------------------
    /* Neumann OKP - flow - id: settings > neumann flow_per_square */
    {
        vector<double>* neumannFlowPS = new vector<double>();
        string tmpPath[] = {"settings", "neumann", "flow_per_squre"};
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), neumannFlowPS)) {
            cout << "   settings: neumann: flow_per_square [] are not defined" << endl;
        }

        for (vector<double>::const_iterator iterNeumann = neumannFlowPS->begin(); iterNeumann < neumannFlowPS->end(); iterNeumann += 2) {
            int tag = (int) (*iterNeumann);
            double value = (double) (*(iterNeumann + 1));

            list<BoundaryData*>::iterator iter = listBoundaryData->end();
            int id = listBoundaryData->size();
            int type = NEUMANN_PER_SQUARE;

            listBoundaryData->insert(iter, new BoundaryData(id, type, value, tag));
        }
    }
    // -------------------------------------------------------------------------
    /* Neumann OKP - flow - id: settings > neumann flow_sum */
    {
        vector<double>* neumannFlowSum = new vector<double>();
        string tmpPath[] = {"settings", "neumann", "flow_sum"};
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), neumannFlowSum)) {
            cout << "   settings: neumann: flow_sum [] are not defined" << endl;
        }

        for (vector<double>::const_iterator iterNeumann = neumannFlowSum->begin(); iterNeumann < neumannFlowSum->end(); iterNeumann += 2) {
            int tag = (int) (*iterNeumann);
            double value = (double) (*(iterNeumann + 1));

            list<BoundaryData*>::iterator iter = listBoundaryData->end();
            int id = listBoundaryData->size();
            int type = NEUMANN_SUM;

            listBoundaryData->insert(iter, new BoundaryData(id, type, value, tag));
        }
    }
    // -------------------------------------------------------------------------
    /* Neumann z vysledku OKP - id: settings > neumannResults */
    {
        vector<int>* neumannResSurfaces = new vector<int>();
        string tmpPath[] = {"settings", "neumann", "Results"};
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), neumannResSurfaces)) {
            cout << "   settings: neumann: Results: []  are not defined" << endl;
        } else {
            cout << "This Boundary Condition is NOT implemented yet." << endl;
        }
    }
    // -------------------------------------------------------------------------
    /* Newton OKP - piezometric head - id: settings > newtonPiezometricHead */
    {
        vector<double>* newtonPiezometricHead = new vector<double>();
        string tmpPath[] = {"settings", "newton", "piezometricHeadW"};
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), newtonPiezometricHead)) {
            cout << "   settings: newton: piezometricHead [] are not defined" << endl;
        }

        for (vector<double>::const_iterator iterNewton = newtonPiezometricHead->begin(); iterNewton < newtonPiezometricHead->end(); iterNewton += 3) {
            int tag = (int) (*iterNewton);
            double value = (double) (*(iterNewton + 1));
            double sigma = (double) (*(iterNewton + 2));

            list<BoundaryData*>::iterator iter = listBoundaryData->end();
            int id = listBoundaryData->size();
            int type = NEWTON_PIEZOMETRIC_HEAD;

            listBoundaryData->insert(iter, new BoundaryData(id, type, value, sigma, tag));
        }
    }
    // -------------------------------------------------------------------------
    /* Newton OKP - headPressure - id: settings > newtonHeadPressure */
    {
        vector<double>* newtonHeadPressure = new vector<double>();
        string tmpPath[] = {"settings", "newton", "headPressureW"};
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), newtonHeadPressure)) {
            cout << "   settings: newton: headPressure [] are not defined" << endl;
        }

        for (vector<double>::const_iterator iterNewton = newtonHeadPressure->begin(); iterNewton < newtonHeadPressure->end(); iterNewton += 3) {
            int tag = (int) (*iterNewton);
            double value = (double) (*(iterNewton + 1));
            double sigma = (double) (*(iterNewton + 2));

            list<BoundaryData*>::iterator iter = listBoundaryData->end();
            int id = listBoundaryData->size();
            int type = NEWTON_HEAD_PRESSURE;

            listBoundaryData->insert(iter, new BoundaryData(id, type, value, sigma, tag));
        }
    }
    // =========================================================================
    /* jmeno vystupniho souboru s Flow ok.p. - id: output > fbc > path */
    string outputFBc;
    {
        string tmpPath[] = {"output", "fbc", "path"};
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), &outputFBc)) {
            cout << "   output: fbc <fileName> is not defined" << endl;
        }
    }
    // -------------------------------------------------------------------------
    /* jmeno vystupniho souboru s Transportni ok.p. - id: output > tbc > path */
    string outputTBc;
    {
        string tmpPath[] = {"output", "tbc", "path"};
        if (!yamlReader.getValue(vector<string > (tmpPath, tmpPath + sizeof (tmpPath) / sizeof (string)), &outputTBc)) {
            cout << "   output: tbc <fileName> is not defined" << endl;
        }
    }
    // =========================================================================
    cout << "Read YML file - End" << endl;
    // =========================================================================
    // =====                                                               =====
    // =========================================================================



    // =========================================================================
    //   read TERRAIN
    // =========================================================================
    list< vector<double>* >* listTerrain = NULL;
    // -------------------------------------------------------------------------
    if (terrainName.length() != 0) {
        listTerrain = new list< vector<double>* >();
        readTerrain(terrainName, listTerrain);
    }

    // =========================================================================
    //   read MESH
    // =========================================================================
    char* meshFormat = (char*) malloc(sizeof (char) * (256));
    list<Physical*>* listPhysical = new list<Physical*>();
    bool* bChangePhysicalId = (bool*) malloc(sizeof (bool));
    *bChangePhysicalId = false;
    // -------------------------------------------------------------------------
    map<int, Node*>* mapNode = new map<int, Node*>();
    map<int, Element*>* mapElement = new map<int, Element*>();
    int numExcludedElm = 0;
    // -------------------------------------------------------------------------
    readMesh(inputMesh, meshFormat, listPhysical, bChangePhysicalId, mapNode, mapElement);

    // =========================================================================
    //   read MATERIAL
    // =========================================================================
    map<int, map<char*, double>* >* mapMaterial = NULL;
    // -------------------------------------------------------------------------
    if (inputMaterial.length() != 0) {
        mapMaterial = new map<int, map<char*, double>* >();
        readMaterial(inputMaterial, mapMaterial);
    }

    // =========================================================================
    //   debug write Elements
    // =========================================================================
    if (false) {
        debugWriteMesh(mapElement);
    }

    // =========================================================================
    //   node modifications by terraint 
    // =========================================================================
    if (terrainName.length() != 0) {
        terrainNodeModification(mapNode, listTerrain);
    }

    // =========================================================================
    //   count Element mass point
    // =========================================================================
    countElementMassPoint(mapElement);

    // =========================================================================
    //   mark fictive elemetns
    // =========================================================================
    markFictiveElements(mapElement, PHYSICAL_GROUP, excludePhysicalIds, &numExcludedElm);
    markFictiveElements(mapElement, ELEMENT_GROUP, excludeElementaryIds, &numExcludedElm);
    // =========================================================================

    // =========================================================================
    //   write MESH
    // =========================================================================
    writeMesh(outputMesh, meshFormat, listPhysical, mapNode, mapElement, excludePhysicalIds, numExcludedElm);

    // =========================================================================
    //   write TIC
    // =========================================================================
    if (outputTic.length() != 0) {
        writeTic(outputTic, mapElement, excludePhysicalIds, numExcludedElm);
    }

    // =========================================================================
    //   read NGH
    // =========================================================================
    map<int, list< int* >* >* mapElmNgh = new map<int, list< int* >* >();
    // -------------------------------------------------------------------------
    readNgh(inputNgh, mapElmNgh);

    // =========================================================================
    //   debug write NGH
    // =========================================================================
    if (false) {
        debugWriteNgh(mapElmNgh);
    }

    // =========================================================================
    //   create Groups of Elemets - PHYSICAL
    // =========================================================================
    map<int, list<int>* >* mapGroupPHYS = new map<int, list<int>* >();
    createGroupElements(mapElement, PHYSICAL_GROUP, mapGroupPHYS);
    writeGroupElements(mapGroupPHYS, mapElement);

    // =========================================================================
    //   create Groups of Elemets - ELEMENT
    // =========================================================================
    map<int, list<int>* >* mapGroupELEM = new map<int, list<int>* >();
    createGroupElements(mapElement, ELEMENT_GROUP, mapGroupELEM);
    writeGroupElements(mapGroupELEM, mapElement);

    // =========================================================================
    //   Boundary Condition calculation
    // =========================================================================
    list<Boundary*>* listBoundary = new list<Boundary*>();
    // -------------------------------------------------------------------------
    createBoundaryConditions(listBoundaryData, mapGroupPHYS, mapElement, mapMaterial, mapElmNgh, listBoundary);


    // =========================================================================
    //   write FBC
    // =========================================================================
    if (outputFBc.length() != 0) {
        writeFBc(outputFBc, listBoundary);
    }

    // =========================================================================
    // Smazani vsech neumannu na elementech, ktere nakonec v siti nebudou - Begin
    // =========================================================================
    /*
    int bcdNeumannSize2 = 10;
    int bcdNeumannCounter2 = 0;
    vector<Neumann> neumannV2(bcdNeumannSize2, Neumann());

    for (map<int, Element*>::const_iterator itElement = mapElement->begin(); itElement != mapElement->end(); itElement++) {
        Element* element = itElement->second;

        for (vector<int>::const_iterator exPhysIds = excludePhysicalIds->begin(); exPhysIds < excludePhysicalIds->end(); exPhysIds++) {
            if (element->getTag(0) == (*exPhysIds)) {
                for (vector<Neumann>::const_iterator neumannPryc = neumannBC.begin(); neumannPryc < neumannBC.end(); neumannPryc++) {
                    if ((*neumannPryc).elmId == element->getLabel()) {
                        // tento pryc
                        (*neumannPryc).elmId = -1;
                    }
                }
            }
        }
    }

    for (vector<Neumann>::const_iterator neumannPryc = neumannBC.begin(); neumannPryc < neumannBC.end(); neumannPryc++) {
        if ((*neumannPryc).elmId != -1) {
            if (bcdNeumannCounter2 == bcdNeumannSize2) {
                bcdNeumannSize2 += 10;
                neumannV2.resize(bcdNeumannSize2, Neumann());
            }
            //    cout << bcdNeumannCounter2 << endl;
            neumannV2[bcdNeumannCounter2++] = (*neumannPryc);
        }
    }
    neumannV2.resize(bcdNeumannCounter2);

    cout << "Pocet vyrazenych neumannu: " << bcdNeumannCounter2 << endl;
     */
    // =========================================================================
    // Smazani vsech neumannu na elementech, ktere nakonec v siti nebudou - End
    // =========================================================================

    // =========================================================================
    //   read Results
    // =========================================================================
    if (inputResults.length() != 0) {
        readResults(inputResults);
    }

    // =========================================================================
    //   write TBC
    // =========================================================================
    if (outputTBc.length() != 0) {
        writeTBc(outputTBc);
    }

    /**************************************************************************/
    /***   Clean MAIN data structures - begin                                 */
    /**************************************************************************/
    // uvolneni alokovane pameti
    free(meshFormat);
    // -------------------------------------------------------------------------
    free(bChangePhysicalId);
    // =========================================================================
    // Smazat obsah a potom i kontejner
    // -------------------------------------------------------------------------
    if (listTerrain != NULL) {
        for (list<vector<double>*>::iterator itTerrain = listTerrain->begin(); itTerrain != listTerrain->end(); itTerrain = listTerrain->begin()) {
            vector<double>* vectorTmp = *itTerrain;
            listTerrain->erase(itTerrain);
            delete vectorTmp;
        }
        delete listTerrain;
    }
    // -------------------------------------------------------------------------
    for (map<int, Node*>::const_iterator itNode = mapNode->begin(); itNode != mapNode->end(); itNode++) {
        Node* node = itNode->second;
        mapNode->erase(itNode->first);
        delete node;
    }
    delete mapNode;
    // -------------------------------------------------------------------------
    for (map<int, Element*>::const_iterator itElement = mapElement->begin(); itElement != mapElement->end(); itElement++) {
        Element* element = itElement->second;
        mapElement->erase(itElement->first);
        delete element;
    }
    delete mapElement;
    // -------------------------------------------------------------------------
    for (list<Physical*>::iterator itPhysical = listPhysical->begin(); itPhysical != listPhysical->end(); itPhysical = listPhysical->begin()) {
        Physical* physical = *itPhysical;
        listPhysical->erase(itPhysical);
        delete physical;
    }
    delete listPhysical;
    // -------------------------------------------------------------------------
    if (mapMaterial != NULL) {
        for (map<int, map<char*, double>* >::const_iterator itSubMap = mapMaterial->begin(); itSubMap != mapMaterial->end(); itSubMap++) {
            map<char*, double>* subMap = itSubMap->second;
            subMap->clear();
            mapMaterial->erase(itSubMap->first);
            delete subMap;
        }
        delete mapMaterial;
    }
    // -------------------------------------------------------------------------
    for (map<int, list< int* >* >::const_iterator itElmNgh = mapElmNgh->begin(); itElmNgh != mapElmNgh->end(); itElmNgh++) {
        list< int* >* listElm = itElmNgh->second;
        for (list<int*>::iterator itList = listElm->begin(); itList != listElm->end(); itList = listElm->begin()) {
            int* arrInt = *itList;
            listElm->erase(itList);
            free(arrInt);
        }
        mapElmNgh->erase(itElmNgh->first);
        delete listElm;
    }
    delete mapElmNgh;
    // -------------------------------------------------------------------------
    if (listBoundary != NULL) {
        for (list<Boundary*>::iterator itBoundary = listBoundary->begin(); itBoundary != listBoundary->end(); itBoundary = listBoundary->begin()) {
            Boundary* boundary = *itBoundary;
            listBoundary->erase(itBoundary);
            delete boundary;
        }
        delete listBoundary;
    }
    /**************************************************************************/
    /***   Clean MAIN data structures - end                                   */
    /**************************************************************************/

    programExit();

    return 0;
}
