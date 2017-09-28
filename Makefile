# gcc.
CXX=	g++

# compiler flags.
CXXFLAGS=	-I/usr/local/include/eigen3
#CXXFLAGS+=	-Ofast
CXXFLAGS+=	-O2 -g2
CXXFLAGS+=	-std=c++11
LDFLAGS=	-lstdc++

CLEANFILES= tools

all:	tools

clean:
	@rm -rf ${CLEANFILES}

