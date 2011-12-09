/* 
 * File:   TElementFactory.h
 * Author: dalibor
 *
 * Created on June 1, 2010, 2:14 PM
 */

#ifndef _TELEMENTFACTORY_H
#define	_TELEMENTFACTORY_H

#include "element.h"
#include "node.h"

class TElementFactory {
public:
    TElement* getElement(int, int, int, int*, TNode**);

    static TElementFactory* getInstance() {
        if (TElementFactory::instance == NULL) {
            TElementFactory::instance = new TElementFactory();
        }
        return TElementFactory::instance;
    }

private:
    TElementFactory();
    ~TElementFactory();

    static TElementFactory* instance;
};

#endif	/* _TELEMENTFACTORY_H */

