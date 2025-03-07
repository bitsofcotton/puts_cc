#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cctype>
#include <random>
#include <assert.h>

#define int int64_t
#include "lieonn.hh"
typedef myfloat num_t;
#include "corpus.hh"
std::vector<std::string> words;

std::vector<std::string> delimiter;
std::vector<std::string> csvelim;
std::vector<std::string> csvdelim;

void usage() {
  std::cout << "tools (lword|lbalance|toc|lack|redig|stat|findroot|diff|same|prep|pred)" << std::endl;
}

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

inline std::string utf8align(const std::string& tob) {
  int head = 0;
  while(head < tob.size() && (tob[head] & 0xc0) == 0x80) head ++;
  int tail = head;
  for(int j = head; j < tob.size(); j ++) if((tob[j] & 0xc0) != 0x80) tail = j;
  if(-- tail <= head) return tob.substr(0, 0);
  int cnt(0);
  for(int j = head; j <= tail; j ++) if((tob[j] & 0xc0) != 0x80) cnt ++;
  return cnt <= 1 ? tob.substr(0, 0) : tob.substr(head, tail - head + 1);
}

inline void makelword(vector<string>& words, const std::string& input, const bool& show = false, const bool& utf8 = true) {
  words.insert(words.end(), csvelim.begin(),  csvelim.end());
  words.insert(words.end(), csvdelim.begin(), csvdelim.end());
  std::sort(words.begin(), words.end());
  words.erase(std::unique(words.begin(), words.end()), words.end());
  std::vector<gram_t<std::string> > found;
  const auto lwords(lword<char, std::string>(int(log(num_t(int(input.size() ))) / log(num_t(int(2)) ) )).compute(input));
  for(auto itr = lwords.begin(); itr != lwords.end(); ++ itr) {
    if(itr->rptr.size() < 2 && itr->str.size() < 3)
      continue;
    const auto lb(std::lower_bound(found.begin(), found.end(), *itr));
    if(found.begin() <= lb && lb < found.end() && lb->str == itr->str)
      lb->rptr.insert(lb->rptr.end(), itr->rptr.begin(), itr->rptr.end());
    else
      found.emplace_back(*itr);
  }
  for(auto itr = found.begin(); itr != found.end(); ++ itr) {   
    std::sort(itr->rptr.begin(), itr->rptr.end());
    itr->rptr.erase(std::unique(itr->rptr.begin(), itr->rptr.end()), itr->rptr.end());
  }
  std::sort(found.begin(), found.end(), lessCount<std::string>);
  found.erase(std::unique(found.begin(), found.end()), found.end());
  words.reserve(words.size() + found.size());
  for(auto itr(found.begin()); itr < found.end(); ++ itr) {
    const auto tob(utf8 ? utf8align(itr->str) : itr->str);
    if(! tob.size()) continue;
    words.emplace_back(tob);
    if(show) std::cout << tob << ", " << itr->rptr.size() << std::endl;
  }
  std::sort(words.begin(), words.end());
  words.erase(std::unique(words.begin(), words.end()), words.end());
  auto mydelim(delimiter);
  mydelim.insert(mydelim.end(), words.begin(), words.end());
  sort(mydelim.begin(), mydelim.end());
  auto inputs(cutText(input, words, mydelim));
  std::sort(inputs.begin(), inputs.end());
  inputs.erase(std::unique(inputs.begin(), inputs.end()), inputs.end());
  if(utf8)
    for(int i = 0; i < inputs.size(); i ++) {
      inputs[i] = utf8align(inputs[i]);
      if(inputs[i].size()) words.emplace_back(inputs[i]);
    }
  else
    words.insert(words.end(), inputs.begin(), inputs.end());
  std::sort(words.begin(), words.end());
  words.erase(std::unique(words.begin(), words.end()), words.end());
  if(show)
    for(int i = 0; i < inputs.size(); i ++)
      if(inputs[i].size()) std::cout << inputs[i] << ", 1" << std::endl;
  return;
}

#undef int
int main(int argc, const char* argv[]) {
#if defined(NOARCFOUR)
  srandom_dev();
#endif
#define int int64_t
  if(argc < 2) {
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
  std::sort(delimiter.begin(), delimiter.end());
  std::sort(csvelim.begin(), csvelim.end());
  std::sort(csvdelim.begin(), csvdelim.end());
  std::string input, line;
  while(std::getline(std::cin, line)) {
    for(int i = 0; i < line.size(); i ++)
      if(line[i] == '<' || line[i] == '>' || line[i] == '&') line[i] = '!';
    input += line + std::string("\n");
  }
  if(strcmp(argv[1], "lword") != 0) {
    if(2 < argc) {
      words = cutText(loadbuf(argv[2]).second, csvelim, csvdelim, true);
      if(! words.size()) makelword(words, input);
    } else
      makelword(words, input);
  } else if(2 < argc)
    words = cutText(loadbuf(argv[2]).second, csvelim, csvdelim, true);
  std::sort(words.begin(), words.end());
  words.erase(std::unique(words.begin(), words.end()), words.end());
  // setup as default parameters.
  const int szwindow(sqrt(num_t(int(input.size()))));
  const int outblock(sqrt(sqrt(num_t(int(input.size() )) )) );
  const num_t redig(int(1));
  // N.B. cbrt is optimal but we use square except for predTOC.
  const int nrwords(pow(num_t(int(19683)), num_t(int(2)) / num_t(int(3)) ));
  if(std::strcmp(argv[1], "lword") == 0)
    makelword(words, input, true);
  else if(std::strcmp(argv[1], "lbalance") == 0) {
    const auto cinput(cutText(input, csvelim, delimiter));
    const auto idxs(pseudoWordsBalance<double, std::string>(cinput, words, szwindow));
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
    preparedTOC<num_t, std::string>(std::cout, input, detailwords, details, tocwords, tocs, delimiter, szwindow, outblock, nrwords, redig, std::strcmp(argv[1], "lack") == 0);
    std::cout << std::endl << "<br/></body></html>";
  } else if(std::strcmp(argv[1], "reconstruct") == 0)
    std::cout << corpus<double, std::string>(input, delimiter).serialize() << std::endl;
  else if(std::strcmp(argv[1], "redig") == 0) {
    std::vector<double> emph;
    emph.push_back(4.);
    emph.push_back(1.);
    emph.push_back(.25);
    for(int ei = 0; ei < emph.size(); ei ++) {
      for(int i = 0; i < input.size() / szwindow; i ++)
        std::cout << corpus<double, std::string>(input.substr(i * szwindow, szwindow), delimiter).reDig(emph[ei]).serialize() << std::endl;
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
    diff<num_t, std::string>(std::cout, input, details, detailwords, details2, detailwords2, delimiter, szwindow, outblock, nrwords, redig, strcmp(argv[1], "same") == 0);
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
    optimizeTOC<num_t, std::string>(std::cout, input, rdetails, rdetailwords, delimiter, szwindow, outblock, nrwords, redig, std::strcmp(argv[1], "findroot") == 0);
    std::cout << "<hr/>" << std::endl << "</body></html>" << std::endl;
  } else if(std::strcmp(argv[1], "pred") == 0) {
    std::vector<std::string> details;
    std::vector<std::string> detailwords;
    for(int iidx = 3; iidx < argc; iidx ++) {
      const auto work(loadbuf(argv[iidx]));
      details.push_back(work.second);
      detailwords.push_back(work.first);
    }
    words.insert(words.end(), detailwords.begin(), detailwords.end());
    std::sort(words.begin(), words.end());
    words.erase(std::unique(words.begin(), words.end()), words.end());
    predTOC<num_t, std::string>(std::cout, input, detailwords, details, delimiter, szwindow, nrwords, redig);
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

