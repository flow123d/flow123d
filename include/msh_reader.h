/* 
 * File:   mesh.h
 * Author: dalibor
 *
 * Created on October 3, 2010, 11:23 AM
 */

#ifndef _MESHREADER_H
#define	_MESHREADER_H

#include "mesh.h"

/**
 * MeshReader is ONLY basic for family of mesh readers.
 * ALL methods MUST BE virtual.
 */
class MeshReader {
private:

public:

    MeshReader() {
    }

    ~MeshReader() {
    }

    /**
     *  Read mesh from file
     */
    virtual void read(const char*, Mesh*) {
        ASSERT(false, "MeshRedaer->read(const char*, Mesh*) is ONLY virtual method.\n");
    }

};

#endif	/* _MESHREADER_H */
