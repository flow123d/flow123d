/* 
 * File:   gmshmeshreader.h
 * Author: dalibor
 *
 * Created on October 3, 2010, 11:23 AM
 */

#ifndef _GMSHMESHREADER_H
#define	_GMSHMESHREADER_H

#include "msh_reader.h"
#include "mesh.h"

class GmshMeshReader : public MeshReader {
private:
    void read_nodes(FILE*, Mesh*);
    void read_elements(FILE*, Mesh*);

    char supported_element_type(int);
    void parse_element_line(ElementVector&, char*, Mesh* mesh);
    void element_type_specific(ElementFullIter);
    void element_allocation_independent(ElementFullIter);

public:
    GmshMeshReader();
    ~GmshMeshReader();

    /**
     *  Read mesh from file
     */
    void read(const char*, Mesh*);

};

#endif	/* _GMSHMESHREADER_H */

