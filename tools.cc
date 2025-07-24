#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cctype>
#include <assert.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#if !defined(_OLDCPP_)
#define int int64_t
#endif
#include "lieonn.hh"
typedef myfloat num_t;
std::vector<std::string> words;

std::pair<std::string, std::string> loadbuf(const char* filename) {
  std::ifstream input;
  std::string   line;
  std::string   inbuf;
  input.open(filename);
  while(getline(input, line)) {
    for(int i = 0; i < line.size(); i ++)
      if(line[i] == '<' || line[i] == '>' || line[i] == '&') line[i] = '!';
    inbuf += line + std::string("\n");
    if(input.eof() || input.bad())
      break;
  }
  input.close();
  
  std::string name0(filename);
  int slash = - 1;
  for(int j = 0; j < name0.size(); j ++)
    if(name0[j] == '/')
      slash = j;
  std::string name;
  slash ++;
  for(int j = slash; j < name0.size(); j ++)
    name += name0[j];
  
  return std::make_pair(name, inbuf);
}

#if !defined(_OLDCPP_)
#undef int
#endif
int main(int argc, const char* argv[]) {
#if !defined(_OLDCPP_)
#define int int64_t
#endif
  std::string input, line;
  std::vector<std::string> delimiter;
  std::vector<std::string> csvelim;
  std::vector<std::string> csvdelim;
  if(argc < 2) goto usage;
  delimiter.push_back(string("."));
  delimiter.push_back(string(","));
  delimiter.push_back(string("\'"));
  delimiter.push_back(string("\""));
  delimiter.push_back(string("。"));
  delimiter.push_back(string("、"));
  delimiter.push_back(string("「"));
  delimiter.push_back(string("」"));
  delimiter.push_back(string("("));
  delimiter.push_back(string(")"));
  csvelim.push_back(string(" "));
  csvelim.push_back(string("\t"));
  csvdelim.push_back(string(","));
  csvdelim.push_back(string("\r"));
  csvdelim.push_back(string("\n"));
  std::sort(delimiter.begin(), delimiter.end());
  std::sort(csvelim.begin(), csvelim.end());
  std::sort(csvdelim.begin(), csvdelim.end());
  while(std::getline(std::cin, line)) {
    for(int i = 0; i < line.size(); i ++)
      if(line[i] == '<' || line[i] == '>' || line[i] == '&') line[i] = '!';
    input += line + std::string("\n");
  }
  words.insert(words.end(), csvelim.begin(),  csvelim.end());
  words.insert(words.end(), csvdelim.begin(), csvdelim.end());
  std::sort(words.begin(), words.end());
  words.erase(std::unique(words.begin(), words.end()), words.end());
  if(strcmp(argv[1], "lword") != 0) {
    if(2 < argc) {
      words = cutText(loadbuf(argv[2]).second, csvelim, csvdelim, true);
      if(! words.size()) makelword<num_t, std::string>(words, input, delimiter);
    } else
      makelword<num_t, std::string>(words, input, delimiter);
  } else if(2 < argc)
    words = cutText(loadbuf(argv[2]).second, csvelim, csvdelim, true);
  std::sort(words.begin(), words.end());
  words.erase(std::unique(words.begin(), words.end()), words.end());
  if(std::strcmp(argv[1], "lword") == 0)
    makelword<num_t, std::string>(words, input, delimiter, true, true, - 1);
  else if(std::strcmp(argv[1], "lbalance") == 0) {
    const vector<string> cinput(cutText(input, csvelim, delimiter));
    const vector<int> idxs(pseudoWordsBalance<num_t, std::string>(cinput, words));
    std::cout << idxs.size() << "sets." << std::endl;
    for(int i = 0; i < idxs.size(); i ++)
      std::cout << cinput[idxs[i]] << std::endl;
  } else if(std::strcmp(argv[1], "toc") == 0 ||
            std::strcmp(argv[1], "lack") == 0) {
    std::vector<std::string> details;
    std::vector<std::string> tocs;
    std::vector<std::string> detailwords;
    std::vector<std::string> tocwords;
    bool toc(false);
    for(int iidx = 3; iidx < argc; iidx ++) {
      if(std::string(argv[iidx]) == std::string("-toc")) {
        toc = true;
        continue;
      }
      const pair<string, string> work(loadbuf(argv[iidx]));
      if(toc) {
        tocs.push_back(work.second);
        tocwords.push_back(work.first);
      } else {
        details.push_back(work.second);
        detailwords.push_back(work.first);
      }
    }
    words.insert(words.end(), tocwords.begin(),    tocwords.end());
    words.insert(words.end(), detailwords.begin(), detailwords.end());
    std::sort(words.begin(), words.end());
    words.erase(std::unique(words.begin(), words.end()), words.end());
    std::cout << "<html><head><link rel=\"stylesheet\" type=\"text/css\" href=\"../../style.css\"><meta charset=\"utf-8\" /></head>" << std::endl;
    std::cout << "<body>";
    preparedTOC<num_t, std::string>(std::cout, input, detailwords, details, tocwords, tocs, delimiter, std::strcmp(argv[1], "lack") == 0);
    std::cout << std::endl << "<br/></body></html>";
  } else if(std::strcmp(argv[1], "reconstruct") == 0)
    std::cout << corpus<num_t, std::string>(input, delimiter).serialize() << std::endl;
  else if(std::strcmp(argv[1], "redig") == 0) {
    std::vector<num_t> emph;
    emph.push_back(4.);
    emph.push_back(1.);
    emph.push_back(.25);
    const int szwindow(sqrt(num_t(int(input.size()))));
    for(int ei = 0; ei < emph.size(); ei ++) {
      for(int i = 0; i < input.size() / szwindow; i ++)
        std::cout << corpus<num_t, std::string>(input.substr(i * szwindow, szwindow), delimiter).reDig(emph[ei]).serialize() << std::endl;
      std::cout << std::endl << std::endl;
    }
  } else if(std::strcmp(argv[1], "diff") == 0 ||
            std::strcmp(argv[1], "same") == 0) {
    std::vector<std::string> details, details2;
    std::vector<std::string> detailwords, detailwords2;
    bool second(false);
    for(int iidx = 3; iidx < argc; iidx ++) {
      if(std::string(argv[iidx]) == std::string("-dict")) {
        second = false;
        continue;
      } else if(std::string(argv[iidx]) == std::string("-dict2")) {
        second = true;
        continue;
      }
      const pair<string, string> work(loadbuf(argv[iidx]));
      if(second) {
        details2.push_back(work.second);
        detailwords2.push_back(work.first);
      } else {
        details.push_back(work.second);
        detailwords.push_back(work.first);
      }
    }
    words.insert(words.end(), detailwords.begin(),  detailwords.end());
    words.insert(words.end(), detailwords2.begin(), detailwords2.end());
    std::sort(words.begin(), words.end());
    words.erase(std::unique(words.begin(), words.end()), words.end());
    std::cout << "<html><head><link rel=\"stylesheet\" type=\"text/css\" href=\"../../style.css\"><meta charset=\"utf-8\" /></head>" << std::endl;
    std::cout << "<body>";
    diff<num_t, std::string>(std::cout, input, details, detailwords, details2, detailwords2, delimiter, strcmp(argv[1], "same") == 0);
    std::cout << "<hr/>" << std::endl << "</body></html>" << std::endl;
  } else if(std::strcmp(argv[1], "stat") == 0 ||
            std::strcmp(argv[1], "findroot") == 0) {
    std::vector<std::string> rdetails;
    std::vector<std::string> rdetailwords;
    for(int iidx = 3; iidx < argc; iidx ++) {
      const pair<string, string> work(loadbuf(argv[iidx]));
      rdetails.push_back(work.second);
      rdetailwords.push_back(work.first);
    }
    words.insert(words.end(), rdetailwords.begin(), rdetailwords.end());
    std::sort(words.begin(), words.end());
    words.erase(std::unique(words.begin(), words.end()), words.end());
    std::cout << "<html><head><link rel=\"stylesheet\" type=\"text/css\" href=\"../../style.css\"><meta charset=\"utf-8\" /></head>" << std::endl;
    std::cout << "<body>";
    optimizeTOC<num_t, std::string>(std::cout, input, rdetails, rdetailwords, delimiter, std::strcmp(argv[1], "findroot") == 0);
    std::cout << "<hr/>" << std::endl << "</body></html>" << std::endl;
  } else if(std::strcmp(argv[1], "pred") == 0) {
    std::vector<std::string> details;
    std::vector<std::string> detailwords;
    for(int iidx = 3; iidx < argc; iidx ++) {
      const pair<string, string> work(loadbuf(argv[iidx]));
      details.push_back(work.second);
      detailwords.push_back(work.first);
    }
    words.insert(words.end(), detailwords.begin(), detailwords.end());
    std::sort(words.begin(), words.end());
    words.erase(std::unique(words.begin(), words.end()), words.end());
    predTOC<num_t, std::string>(std::cout, input, detailwords, details, delimiter);
  } else if(std::strcmp(argv[1], "prep") == 0) {
    std::vector<std::string> buf;
    const int szwindow(sqrt(num_t(int(input.size()))));
    for(int i = 0; i < input.size() / szwindow + 1; i ++) {
      const vector<string> work(corpus<num_t, std::string>(input.substr(i * szwindow, std::min(szwindow, int(input.size()) - i * szwindow)), delimiter).reverseLink().first);
      buf.insert(buf.end(), work.begin(), work.end());
    }
    std::sort(buf.begin(), buf.end());
    buf.erase(std::unique(buf.begin(), buf.end()), buf.end());
    for(int i = 0; i < buf.size(); i ++)
      std::cout << buf[i] << std::endl;
  } else if(std::strcmp(argv[1], "optdict") == 0)
    assert(0 && "group dicts: not implemented around NOT word tables.");
  else if(std::strcmp(argv[1], "conflict") == 0)
    assert(0 && "conflict : not implemented around NOT word tables.");
  else if(std::strcmp(argv[1], "negate") == 0)
    assert(0 && "negate: not implemented around NOT word tables.");
  else if(std::strcmp(argv[1], "consistency") == 0)
    assert(0 && "consistancy : Logics check so far...");
  else if(std::strcmp(argv[1], "logiccheck") == 0)
    assert(0 && "logic check : Logics check so far...");
  else goto usage;
  return 0;
 usage:
  std::cout << "tools (lword|lbalance|toc|lack|redig|stat|findroot|diff|same|prep|pred)" << std::endl;
  return - 2;
}

