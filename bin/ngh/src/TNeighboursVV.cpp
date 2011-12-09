/* 
 * File:   TNeighboursVV.cpp
 * Author: dalibor
 * 
 * Created on 2. červen 2010, 17:20
 */

#include "TNeighboursVV.h"

#include <iostream>
#include "element.h"
#include "neighbour.h"

TNeighboursVV::TNeighboursVV() {
}

TNeighboursVV::~TNeighboursVV() {
}

void TNeighboursVV::create(TProblem* problem, TMesh* mesh) {
    std::cout << "  Creating VV neighbours... ";

    TElement* ele1;
    TElement* ele2;
    TNeighbour* ngh;

    double coef;
    int count = 0;
    int pc = 0, p = 0;

    int pinc;

    int perc = mesh->getNumElements() / 100;
    if (perc < 100) {
        pinc = 100 / mesh->getNumElements();
    } else {
        pinc = 1;
    }

    FOR_ELEMENTS(ele1, mesh->getFirstElement()) {
        //      std::cout << ele1->GetId() << "\n";
        for (ele2 = ele1->getNext(); ele2 != NULL; ele2 = ele2->getNext()) {
            //      if (ele1->GetId() == 36)
            if (ele1->getLabel() == 36 && ele2->getLabel() == 42) {
                std::cout << ele1->getLabel() << " " << ele2->getLabel() << "\n";
            }

            if (ele1->AreVVNeighbours(problem, ele2, coef)) {
                ngh = new TNeighbour(ele1, ele2, coef);
                mesh->AddNeighbour(ngh);
                count++;
            }
        }
        if (pc == perc || ele1->getNext() == NULL) {
            p += pinc;
            if (ele1->getNext() == NULL) {
                p = 100;
            }
            if (count != 0) {
                std::cout << "\r  Creating VV neighbours... " << p << "% ";
            }
            pc = -1;
        }
        pc++;
    }
    std::cout << count << " neighbours created. OK.\n";
    return;
}
