CXX=	clang++

# compiler flags.
CXXFLAGS+=	-I/usr/local/include/eigen3
CXXFLAGS+=	-Ofast -mtune=native -g2
CXXFLAGS+=	-std=c++11
LDFLAGS+=	-lc++

# without eigen, this lacks comparestructure function.
# CXXFLAGS+=	-D_WITHOUT_EIGEN_

CLEANFILES= *.o tools

all:	tools

clean:
	@rm -rf ${CLEANFILES}

tools.o:        tools.cc sparse.hh corpus.hh file2eigen.hh lword.hh
test.o:		test.cc  sparse.hh corpus.hh file2eigen.hh lword.hh

