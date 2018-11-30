CXX=	clang++

# compiler flags.
CXXFLAGS+=	-I/usr/local/include/eigen3
CXXFLAGS+=	-Ofast -g0
#CXXFLAGS+=	-O2 -g2
CXXFLAGS+=	-std=c++11
LDFLAGS+=	-lc++

# assertion
# CXXFLAGS+=	-D_STRICT_WORD_ASSERT_

CLEANFILES= *.o tools test

all:	tools test

clean:
	@rm -rf ${CLEANFILES}

tools.o:        tools.cc sparse.hh corpus.hh file2eigen.hh lword.hh
test.o:		test.cc  sparse.hh corpus.hh file2eigen.hh lword.hh

