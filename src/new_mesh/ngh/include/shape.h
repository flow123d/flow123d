
#ifndef _SHAPE_H
#define	_SHAPE_H

#ifdef	__cplusplus
extern "C" {
#endif

    typedef enum Shapes {
        UNKNOWN = 1000,
        POINT = 0,
        LINE = 1,
        TRIANGLE = 2,
        TETRAHEDRON = 4
    } TShape;


#ifdef	__cplusplus
}
#endif

#endif	/* _SHAPE_H */
