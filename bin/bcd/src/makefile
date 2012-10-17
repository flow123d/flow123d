#include ../../../makefile.in

CINCLUDE=-I../include

OBJS=\
bcd.o\
condition.o\
mesh.o\
output.o\
system.o\
topology_data.o

%.o: %.c
	${CC} -c ${CFLAGS} ${CINCLUDE} $<

all: ${OBJS}
	${CC} ${CFLAGS} -o ../bin/bcd ${OBJS} -lm 

clean:
	rm -f *.o
	rm -f ../bin/bcd
