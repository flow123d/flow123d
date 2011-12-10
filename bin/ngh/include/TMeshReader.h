/* 
 * File:   TMeshReader.h
 * Author: dalibor
 *
 * Created on June 1, 2010, 12:09 PM
 */

#ifndef _TMESHREADER_H
#define	_TMESHREADER_H

#include "mesh.h"

class TMeshReader {
public:
    TMeshReader();
    ~TMeshReader();

    void read(char*, TMesh*);

private:
    void seekNodes(FILE*);
    TNode* parseNodeLine(char*);
    void readNodes(FILE*, TMesh*);

    void seekElements(FILE*);
    TElement* parseElementLine(char*, TMesh*);
    void readElements(FILE*, TMesh*);
};

#endif	/* _TMESHREADER_H */

