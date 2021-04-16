#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <locale>
#include <algorithm>
#include <vector>
#include <map>
#include <utility>
#include <assert.h>

#include "ifloat.hh"
#define int int64_t
typedef myfloat num_t;
#include "simplelin.hh"
#include "corpus.hh"

std::vector<std::string> words;

void usage() {
  std::cout << "tools (lword|lbalance|toc|lack|redig|stat|findroot|diff|same|prep)" << std::endl;
}

//const int    szwindow(120);
const int    szwindow(1500);
const int    Mbalance(40);
const double scorethresh(.25);
const double dscorethresh(.001);
const double threshin(0.);
const double redig(1.1);
std::vector<std::string> delimiter;
std::vector<std::string> csvelim;
std::vector<std::string> csvdelim;

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

#undef int
int main(int argc, const char* argv[]) {
#define int int64_t
  if(argc < 3) {
    usage();
    return - 2;
  }
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
  words = cutText(loadbuf(argv[2]).second, csvelim, csvdelim, true);
  std::string input, line;
  while(std::getline(std::cin, line)) {
    for(int i = 0; i < line.size(); i ++)
      if(line[i] == '<' || line[i] == '>' || line[i] == '&') line[i] = '!';
    input += line + std::string("\n");
  }
  if(std::strcmp(argv[1], "lword") == 0) {
    words.insert(words.end(), csvelim.begin(),  csvelim.end());
    words.insert(words.end(), csvdelim.begin(), csvdelim.end());
    std::vector<gram_t<std::string> > found;
    for(int i = 2; i < 80; i ++) {
      const auto lwords(lword<char, std::string>(80, i).compute(input));
      for(auto itr = lwords.begin(); itr != lwords.end(); ++ itr) {
        if(itr->rptr.size() < 2 && itr->str.size() < 3)
          continue;
        const auto lb(std::lower_bound(found.begin(), found.end(), *itr, lessCount<std::string>));
        if(found.begin() <= lb && lb < found.end() && lb->str == itr->str)
          lb->rptr.insert(lb->rptr.end(), itr->rptr.begin(), itr->rptr.end());
        else
          found.emplace_back(*itr);
      }
    }
    for(auto itr = found.begin(); itr != found.end(); ++ itr) {   
      std::sort(itr->rptr.begin(), itr->rptr.end());
      itr->rptr.erase(std::unique(itr->rptr.begin(), itr->rptr.end()), itr->rptr.end());
    }
    std::sort(found.begin(), found.end(), lessCount<std::string>);
    found.erase(std::unique(found.begin(), found.end()), found.end());
    words.reserve(words.size() + found.size());
    for(auto itr(found.begin()); itr < found.end(); ++ itr) {
      const auto& tob(itr->str);
      words.emplace_back(tob);
      std::cout << tob << ", " << itr->rptr.size() << std::endl;
    }
    std::sort(words.begin(), words.end());
    words.erase(std::unique(words.begin(), words.end()), words.end());
    auto inputs(cutText(input, words, delimiter));
    std::sort(inputs.begin(), inputs.end());
    inputs.erase(std::unique(inputs.begin(), inputs.end()), inputs.end());
    for(int i = 0; i < inputs.size(); i ++)
      std::cout << inputs[i] << ", 1" << std::endl;
  } else if(std::strcmp(argv[1], "lbalance") == 0) {
    const auto cinput(cutText(input, csvelim, delimiter));
    const auto idxs(pseudoWordsBalance<double, std::string>(cinput, words, Mbalance));
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
      const auto work(loadbuf(argv[iidx]));
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
    preparedTOC<double, std::string>(std::cout, input, detailwords, details, tocwords, tocs, delimiter, szwindow, - scorethresh, threshin, redig, std::strcmp(argv[1], "lack") == 0);
    std::cout << std::endl << "<br/></body></html>";
  } else if(std::strcmp(argv[1], "reconstruct") == 0)
    std::cout << corpus<double, std::string>(input, delimiter).serialize() << std::endl;
  else if(std::strcmp(argv[1], "redig") == 0) {
    std::vector<double> emph;
    emph.push_back(4.);
    emph.push_back(1.);
    emph.push_back(.25);
    for(int ei = 0; ei < emph.size(); ei ++) {
      for(int i = 0; i < input.size() / szwindow + 1; i ++)
        std::cout << corpus<double, std::string>(input.substr(i * szwindow, std::min(szwindow, int(input.size()) - i * szwindow)), delimiter).reDig(emph[ei]).serialize() << std::endl;
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
      const auto work(loadbuf(argv[iidx]));
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
    diff<double, std::string>(std::cout, input, details, detailwords, details2, detailwords2, delimiter, szwindow, - dscorethresh, threshin, redig, strcmp(argv[1], "same") == 0);
    std::cout << "<hr/>" << std::endl << "</body></html>" << std::endl;
  } else if(std::strcmp(argv[1], "stat") == 0 ||
            std::strcmp(argv[1], "findroot") == 0) {
    std::vector<std::string> rdetails;
    std::vector<std::string> rdetailwords;
    for(int iidx = 3; iidx < argc; iidx ++) {
      const auto work(loadbuf(argv[iidx]));
      rdetails.push_back(work.second);
      rdetailwords.push_back(work.first);
    }
    words.insert(words.end(), rdetailwords.begin(), rdetailwords.end());
    std::sort(words.begin(), words.end());
    words.erase(std::unique(words.begin(), words.end()), words.end());
    std::cout << "<html><head><link rel=\"stylesheet\" type=\"text/css\" href=\"../../style.css\"><meta charset=\"utf-8\" /></head>" << std::endl;
    std::cout << "<body>";
    optimizeTOC<double, std::string>(std::cout, input, rdetails, rdetailwords, delimiter, szwindow, - scorethresh, threshin, redig, std::strcmp(argv[1], "findroot") == 0);
    std::cout << "<hr/>" << std::endl << "</body></html>" << std::endl;
  } else if(std::strcmp(argv[1], "prep") == 0) {
    std::vector<std::string> buf;
    for(int i = 0; i < input.size() / szwindow + 1; i ++) {
      const auto work(corpus<double, std::string>(input.substr(i * szwindow, std::min(szwindow, int(input.size()) - i * szwindow)), delimiter).reverseLink().first);
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
  else {
    usage();
    return - 2;
  }
  return 0;
}

