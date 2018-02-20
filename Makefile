CXX=	clang++

# compiler flags.
CXXFLAGS+=	-I/usr/local/include/eigen3
CXXFLAGS+=	-Ofast -g0
#CXXFLAGS+=	-O2 -g2
CXXFLAGS+=	-std=c++11
LDFLAGS+=	-lstdc++

CLEANFILES= *.o tools test

all:	tools test

clean:
	@rm -rf ${CLEANFILES}

tools.o:        tools.cc corpus.hh file2eigen.hh lword.hh
test.o:		test.cc corpus.hh file2eigen.hh lword.hh

