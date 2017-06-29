# gcc.
CXX=	g++

# compiler flags.
CXXFLAGS=	-I/usr/local/include/eigen3
#CXXFLAGS+=	-Ofast
CXXFLAGS+=	-O2 -g2
CXXFLAGS+=	-std=c++14
LDFLAGS=	-lstdc++

CLEANFILES= *.o *.dSYM tools

all:	tools
tools:	corpus.o tools.o

clean:
	@rm -rf ${CLEANFILES}

