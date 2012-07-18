
CINCLUDE=-I../include -I./simpleini

OBJS=\
TNeighboursVV.o\
TElementFactory.o\
TMeshReader.o\
TInstanceProfiler.o\
abscissa.o\
bisector.o\
edge.o\
element.o\
intersection.o\
main.o\
mathfce.o\
matrix.o\
mesh.o\
neighbour.o\
node.o\
output.o\
plain.o\
point.o\
polygon.o\
problem.o\
side.o\
system.o\
tetrahedron.o\
triangle.o\
vector.o\
vertex.o

%.o: %.cpp
	${CXX} -c ${CFLAGS} ${CINCLUDE} $<

all: ${OBJS}
	${CXX} ${CFLAGS} -o ../bin/ngh ${OBJS} -lm 

clean:
	rm -f *.o
	rm -f ../../bin/ngh
