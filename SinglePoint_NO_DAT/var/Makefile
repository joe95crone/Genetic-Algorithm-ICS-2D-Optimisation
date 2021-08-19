OBJECTS=config_parser.o main.o pisa_io.o pop.o problem_mgr.o var_app.o
CFLAGS=-c -ggdb -O0
CC=g++

var: ${OBJECTS}
	${CC} ${OBJECTS} -o $@

.cpp.o:
	${CC} ${CFLAGS} $< -o $@

clean:
	rm -f var *.o *~
