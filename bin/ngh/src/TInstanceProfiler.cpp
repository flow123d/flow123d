
#include <stdio.h>

#include "TInstanceProfiler.h"
#include "side.h"
#include "point.h"
#include "abscissa.h"
#include "matrix.h"
#include "element.h"
#include "myvector.h"
#include "neighbour.h"
#include "node.h"
#include "plain.h"
#include "polygon.h"
#include "tetrahedron.h"
#include "triangle.h"
#include "vertex.h"

void TInstanceProfiler::printResult() {
    char tmp[11];

    std::cout << "***********************************************************\n";
    std::cout << "*****           Instance profiling:                   *****\n";
    std::cout << "*****                                                 *****\n";

    std::cout << "*****     .......................................     *****\n";

    sprintf(tmp, "%10d", TPoint::getNumInstances());
    std::cout << "*****     Created " << tmp << " instances of TPoint      *****\n";
    sprintf(tmp, "%10d", TNode::getNumInstances());
    std::cout << "*****     Created " << tmp << " instances of TNode       *****\n";

    std::cout << "*****     .......................................     *****\n";

    sprintf(tmp, "%10d", TAbscissa::getNumInstances());
    std::cout << "*****     Created " << tmp << " instances of TAbscissa   *****\n";

    sprintf(tmp, "%10d", TBisector::getNumInstances());
    std::cout << "*****     Created " << tmp << " instances of TBisector   *****\n";

    std::cout << "*****     .......................................     *****\n";

    sprintf(tmp, "%10d", TEdge::getNumInstances());
    std::cout << "*****     Created " << tmp << " instances of TEdge       *****\n";

    sprintf(tmp, "%10d", TSide::getNumInstances());
    std::cout << "*****     Created " << tmp << " instances of TSide       *****\n";

    sprintf(tmp, "%10d", TElement::getNumInstances());
    std::cout << "*****     Created " << tmp << " instances of TElement    *****\n";

    sprintf(tmp, "%10d", TVector::getNumInstances());
    std::cout << "*****     Created " << tmp << " instances of TVector     *****\n";

    sprintf(tmp, "%10d", TNeighbour::getNumInstances());
    std::cout << "*****     Created " << tmp << " instances of TNeighbour  *****\n";

    sprintf(tmp, "%10d", TPlain::getNumInstances());
    std::cout << "*****     Created " << tmp << " instances of TPlane      *****\n";

    sprintf(tmp, "%10d", TPolygon::getNumInstances());
    std::cout << "*****     Created " << tmp << " instances of TPolygon    *****\n";

    sprintf(tmp, "%10d", TTetrahedron::getNumInstances());
    std::cout << "*****     Created " << tmp << " instances of TTetrahedron*****\n";

    sprintf(tmp, "%10d", TTriangle::getNumInstances());
    std::cout << "*****     Created " << tmp << " instances of TTriangle   *****\n";

    sprintf(tmp, "%10d", TVertex::getNumInstances());
    std::cout << "*****     Created " << tmp << " instances of TVertex     *****\n";

    std::cout << "***********************************************************\n";
}
