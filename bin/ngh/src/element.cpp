#include "element.h"

#include "config.h"
#include "intersection.h"
#include "mathfce.h"
#include "problem.h"
#include "mesh.h"

using namespace mathfce;

int TElement::numberInstance = 0;

static int NODEARRAY[4][3] = {
    {1, 2, 3},
    {0, 2, 3},
    {0, 1, 3},
    {0, 1, 2}
};

int TElement::generateId() {
    return TElement::numberInstance++;
}

void TElement::Initialize() {
    id = generateId();

    dim = 0;

    n_nodes = 0;
    n_sides = 0;
    n_tags = 0;

    mid = -1;
    rid = -1;

    metrics = 0.0;

    metrics_is_actual = false;

    prev = NULL;
    next = NULL;

    node = NULL;
    side = NULL;

    type = UNKNOWN;

    abscissa = NULL;
    triangle = NULL;
    tetrahedron = NULL;

    return;
}

void TElement::Allocate() {
    switch (type) {
        case LINE:
            n_nodes = 2;
            n_sides = 2;
            break;
        case TRIANGLE:
            n_nodes = 3;
            n_sides = 3;
            break;
        case TETRAHEDRON:
            n_nodes = 4;
            n_sides = 4;
            break;
        default:
            char error[MAXBUFF];
            sprintf(error, "Unknown TYPE (%d) of element (label=%d).", type, label);
            mythrow(error, __LINE__, __FUNC__);
    }

    node = new TNode*[ n_nodes ];
    side = new TSide*[ n_sides ];
}

TElement::TElement(int label, int type, int numTags, int* tags, TNode** nodes) {
    Initialize();

    this->label = label;
    this->type = (TShape) type;

    Allocate();

    n_tags = numTags;
    if (n_tags < 2) {
        char error[MAXBUFF];
        sprintf(error, "Element label=%d, has to have at least two tags.", label);
        mythrow(error, __LINE__, __FUNC__);
    }

    mid = tags[0];
    rid = tags[1];

    int i;

    FOR_ELEMENT_NODES(this, i) {
        this->node[ i ] = nodes[i];
    }
}

TElement::~TElement() {
    if (node != NULL) {
        delete node;
        node = NULL;
    };
    if (side != NULL) {
        delete side;
        side = NULL;
    };
}

int TElement::getLabel() {
    return label;
}

int TElement::GetNNodes() {
    return n_nodes;
}

int TElement::GetNSides() {
    return n_sides;
}

TShape TElement::GetType() {
    return type;
}

int TElement::getNodeLabel(int i) {
    return node[i]->getLabel();
}

TNode* TElement::getNode(int i) {
    return node[i];
}

void TElement::setPrev(TElement* prev) {
    this->prev = prev;
}

void TElement::setNext(TElement* next) {
    this->next = next;
}

TElement* TElement::getNext() {
    return next;
}

void TElement::setSide(int i, TSide* side) {
    this->side[i] = side;
}

void TElement::CreateSides(TMesh* mesh) {
    switch (type) {
        case LINE:
            for (int i = 0; i < 2; ++i) {
                TSide* side = new TSide(this, i, node[ i ]);
                mesh->AddSide(side);
            }
            break;
        case TRIANGLE:
            for (int i = 0; i < 3; ++i) {
                TSide* side = new TSide(this, i, node[ i ], node[ ((i + 1) % 3) ]);
                mesh->AddSide(side);
            }
            break;
        case TETRAHEDRON:
            for (int i = 0; i < 4; ++i) {
                TSide* side = new TSide(this, i, node[NODEARRAY[i][0]], node[ NODEARRAY[i][1]], node[ NODEARRAY[i][2] ]);
                mesh->AddSide(side);
            }
            break;
        default:
            mythrow((char*) "Runtime error - deny point\n", __LINE__, __FUNC__);
    }
}

void TElement::CreateGeometry() {
    switch (type) {
        case LINE:
            abscissa = new TAbscissa(node[0], node[1]);
            break;
        case TRIANGLE:
            triangle = new TTriangle(node[ 0 ], node[ 1 ], node[ 2 ]);
            break;
        case TETRAHEDRON:
            tetrahedron = new TTetrahedron(node[ 0 ], node[ 1 ], node[ 2 ], node[ 3 ]);
            break;
        default:
            mythrow((char*) "Runtime error - deny point\n", __LINE__, __FUNC__);
    }
    return;
}

bool TElement::AreBBNeighbours(TElement* ele2) {
    int si1, si2;

    if (this->type != ele2->type) {
        return false;
    }

    FOR_ELEMENT_SIDES(this, si1) {

        FOR_ELEMENT_SIDES(ele2, si2) {
            if (this->side[ si1 ]->AreBBNeighbours(ele2->side[ si2 ])) {
                return true;
            }
        }
    }
    return false;
}

bool TElement::AreVBNeighbours(TElement* ele2) {
    TElement* ele1 = this;

    if (ele1->type > ele2->type) {
        TElement* etmp = ele1;
        ele1 = ele2;
        ele2 = etmp;
    }

    if (ele1->type == LINE && ele2->type != TRIANGLE) {
        return false;
    }

    if (ele1->type == TRIANGLE && ele2->type != TETRAHEDRON) {
        return false;
    }

    int is;

    FOR_ELEMENT_SIDES(ele2, is) {
        if (ele2->side[ is ]->AreVBNeighbours(ele1)) {
            return true;
        }
    }

    return false;
}

void TElement::Test() {
    char str[MAXBUFF];
    if (IsZero(GetMetrics())) {
        sprintf(str, "Element id = %d is degenerated.", label);
        mythrow(str, __LINE__, __FUNC__);
    }
}

double TElement::GetMetrics() {
    if (metrics_is_actual) {
        return metrics;
    }
    ComputeMetrics();
    return metrics;
}

void TElement::ComputeMetrics() {
    switch (type) {
        case LINE:
            metrics = abscissa->Length();
            break;
        case TRIANGLE:
            metrics = triangle->GetArea();
            break;
        case TETRAHEDRON:
            metrics = tetrahedron->GetVolume();
            break;
        default:
            mythrow((char*) "Runtime error - deny point\n", __LINE__, __FUNC__);
    }
    metrics_is_actual = true;
}

bool TElement::AreVVNeighbours(TProblem* problem, TElement* ele, double& coef) {
    TPoint P;
    TPosition pos;
    TIntersectionType it;

    double t1, t2;

    TElement* ele1;
    TElement* ele2;

    if (type > ele->type) {
        ele1 = ele;
        ele2 = this;
    } else {
        ele1 = this;
        ele2 = ele;
    }

    if (ele1->type == ele2->type && ele1->label > ele2->label) {
        return false;
    }

    switch (ele1->type) {
            //FIRST ELEMENT IS LINE
        case LINE:
            switch (ele2->type) {
                case LINE:
                    if (problem->getCt_11() == ct_ratio1) {
                        if (AreBBNeighbours(ele2)) {
                            return false;
                        }
                        GetIntersection(*ele1->abscissa, *ele2->abscissa, pos, &P);
                        if (pos == intersecting) {
                            coef = 1.0;
                            return true;
                        }
                    }
                    break;
                case TRIANGLE:
                    if (problem->getCt_12() == ct_ratio1) {
                        if (ele1->AreVBNeighbours(ele2)) {
                            return false;
                        }
                        it = unknown;
                        GetIntersection(*ele1->abscissa, *ele2->triangle, it, t1, t2);
                        if (it == point) {
                            coef = 1.0;
                            return true;
                        }
                        if (it == line) {
                            coef = Distance(ele1->abscissa->GetPoint(t1), ele1->abscissa->GetPoint(t2)) /
                                    ele1->abscissa->Length();
                            return true;
                        }
                    }
                    break;
                case TETRAHEDRON:
                    if (problem->getCt_13() == ct_ratio1) {
                        it = unknown;
                        GetIntersection(*ele1->abscissa, *ele2->tetrahedron, it, t1, t2);
                        if (it == line) {
                            coef = Distance(ele1->abscissa->GetPoint(t1), ele1->abscissa->GetPoint(t2)) /
                                    ele1->abscissa->Length();
                            if (IsZero(coef)) return false;
                            return true;
                        }
                    }
                    break;
                default:
                    mythrow((char*) "Runtime error - deny point\n", __LINE__, __FUNC__);
            }
            break;
            //FIRST ELEMENT IS TRIANGLE
        case TRIANGLE:
            switch (ele2->type) {
                case TRIANGLE:
                    if (problem->getCt_22() == ct_ratio1) {
                        if (ele1->AreBBNeighbours(ele2)) {
                            return false;
                        }
                        GetIntersection(*ele1->triangle, *ele2->triangle, it, coef);
                        if ((it == line || it == area) && !IsZero(coef)) {
                            if (it == area) {
                                if (ele1->triangle->GetArea() > ele2->triangle->GetArea()) {
                                    coef = coef / ele1->triangle->GetArea();
                                } else {
                                    coef = coef / ele2->triangle->GetArea();
                                }
                            } else {
                                coef = 1.0;
                            }
                            return true;
                        }
                    }
                    return false;
                case TETRAHEDRON:
                    if (problem->getCt_23() == ct_ratio1) {
                        if (ele1->AreVBNeighbours(ele2)) {
                            return false;
                        }
                        GetIntersection(*ele1->triangle, *ele2->tetrahedron, it, coef);
                        if (it == area && !IsZero(coef)) {
                            coef = coef / ele1->triangle->GetArea();
                            return true;
                        }
                    }
                    return false;
                default:
                    mythrow((char*) "Runtime error - deny point\n", __LINE__, __FUNC__);
            }
        default:
            mythrow((char*) "Runtime error - deny point\n", __LINE__, __FUNC__);
    }
    return false;
}
