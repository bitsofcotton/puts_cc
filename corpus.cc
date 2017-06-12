#include "corpus.hh"

bool partial_compare(const std::string& x0, const std::string& x1) {
  int x0i = 0, x1i = 0;
  char buf = 1;
  while(!buf && x0[x0i] != '\0' && x1[x1i] != '\0') buf = x0[x0i ++] - x1[x1i ++];
  return buf < 0;
}

bool partial_compare_equal(const std::string& x0, const std::string& x1) {
  int x0i = 0, x1i = 0;
  char buf = 0;
  while(!buf && x0[x0i] != '\0' && x1[x1i] != '\0') buf = x0[x0i ++] - x1[x1i ++];
  return !buf;
}

