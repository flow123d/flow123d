/* 
 * File:   TMeshReader.cpp
 * Author: dalibor
 * 
 * Created on June 1, 2010, 12:09 PM
 */

#include "TMeshReader.h"
#include "TElementFactory.h"
#include "system.h"
#include "config.h"

TMeshReader::TMeshReader() {
}

TMeshReader::~TMeshReader() {
}

void TMeshReader::read(char* fileName, TMesh* mesh) {
    FILE* F = fopen(fileName, "rt");
    if (F == NULL) {
        mythrow((char*) "Couldn't open input mesh file.", __LINE__, __FUNC__);
    }

    readNodes(F, mesh);
    readElements(F, mesh);

    fclose(F);
}

void TMeshReader::readNodes(FILE* F, TMesh* mesh) {
    char line[ MAXBUFF ]; // line of data file

    std::cout << "Reading nodes... ";

    seekNodes(F);
    fgets(line, MAXBUFF - 2, F);
    int tmp = atoi(strtok(line, " \t"));
    if (tmp < 1) {
        mythrow((char*) "The number of the nodes has to be greater than 0.", __LINE__, __FUNC__);
    }

    for (int i = 0; i < tmp; i++) {
        fgets(line, MAXBUFF - 2, F);
        TNode* node = parseNodeLine(line);

        mesh->AddNode(node);
    }
    std::cout << tmp << " Nodes read. OK. \n";
}

void TMeshReader::seekNodes(FILE* F) {
    char line[ MAXBUFF ];
    char string[ MAXBUFF ];

    fseek(F, 0, SEEK_SET);

    while (fgets(line, MAXBUFF - 2, F) != NULL) {
        sscanf(line, "%s", string);
        if (strcmpi(string, "$Nodes") == 0) {
            return;
        }
    }
    mythrow((char*) "Cannot find begin of the section of the nodes.", __LINE__, __FUNC__);
}

TNode* TMeshReader::parseNodeLine(char *line) {
    int label = atoi(strtok(line, " \t"));

    double x = atof(strtok(NULL, " \t"));
    double y = atof(strtok(NULL, " \t"));
    double z = atof(strtok(NULL, " \t"));

    TNode* node = new TNode(label, x, y, z);

    return node;
}

void TMeshReader::readElements(FILE* F, TMesh* mesh) {
    char line[ MAXBUFF ]; // line of data file

    std::cout << "Reading elements... ";

    seekElements(F);
    fgets(line, MAXBUFF - 2, F);
    int tmp = atoi(strtok(line, " \t"));
    if (tmp < 1) {
        mythrow((char*) "The number of the elements has to be greater than 0.", __LINE__, __FUNC__);
    }

    for (int i = 0; i < tmp; i++) {
        fgets(line, MAXBUFF - 2, F);
        TElement* element = parseElementLine(line, mesh);

        mesh->AddElement(element);
    }
    std::cout << tmp << " Elements read. OK. \n";
}

void TMeshReader::seekElements(FILE* F) {
    char line[ MAXBUFF ];
    char string[ MAXBUFF ];

    fseek(F, 0, SEEK_SET);

    while (fgets(line, MAXBUFF - 2, F) != NULL) {
        sscanf(line, "%s", string);
        if (strcmpi(string, "$Elements") == 0) {
            return;
        }
    }

    mythrow((char*) "Cannot find begin of the section of the elements.", __LINE__, __FUNC__);
}

TElement* TMeshReader::parseElementLine(char* line, TMesh* mesh) {
    int label = atoi(strtok(line, " \t"));
    int type = atoi(strtok(NULL, " \t"));

    int n_nodes = 0;
    if (type == LINE) {
        n_nodes = 2;
    } else if (type == TRIANGLE) {
        n_nodes = 3;
    } else if (type == TETRAHEDRON) {
        n_nodes = 4;
    }

    int n_tags = atoi(strtok(NULL, " \t"));

    int tmpTagList[n_tags];
    for (int i = 0; i < n_tags; i++) {
        tmpTagList[i] = atoi(strtok(NULL, " \t"));
    }

    TNode * tmpNodeList[4];
    for (int i = 0; i < n_nodes; i++) {
        int nodeLabel = atoi(strtok(NULL, " \t"));
        tmpNodeList[i] = mesh->getNodeLabel(nodeLabel);
    }

    TElement* element = TElementFactory::getInstance()->getElement(type, label, n_tags, tmpTagList, tmpNodeList);

    return element;
}
