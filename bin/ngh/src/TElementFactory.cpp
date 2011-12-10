/* 
 * File:   TElementFactory.cpp
 * Author: dalibor
 * 
 * Created on June 1, 2010, 2:14 PM
 */

#include "TElementFactory.h"

TElementFactory* TElementFactory::instance = NULL;

TElementFactory::TElementFactory() {
}

TElementFactory::~TElementFactory() {
    delete instance;
}

TElement* TElementFactory::getElement(int type, int label, int n_tags, int* tmpTagList, TNode** tmpNodeList) {
    TElement* element = new TElement(label, type, n_tags, tmpTagList, tmpNodeList);

    return element;
}
