#include <Eigen/Core>
#include <Eigen/SVD>
#include <cstdio>
#include <cstring>
#include <string>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <locale>
#include <codecvt>
#include <algorithm>
#include <vector>
#include <map>
#include <utility>

namespace corpus {
#include "lword.hh"
#include "corpus.hh"
#include "file2eigen.hh"
};
using namespace corpus;

void usage() {
  std::cout << "tools (lword|lbalance|corpus|toc|redig|stat|reconstruct|diff)" << std::endl;
}

const int szwindow(200);
const int szblock(8000);
const int Mbalance(40);
const double threshin(.1);
std::vector<std::string> delimiter;
std::vector<std::string> csvelim;
std::vector<std::string> csvdelim;

std::pair<std::string, std::string> loadbuf(const char* filename) {
  std::ifstream input;
  std::string   line;
  std::string   inbuf;
  input.open(filename);
  while(getline(input, line)) {
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

int main(int argc, const char* argv[]) {
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
        auto csv(cutText(loadbuf(argv[3]).second, csvelim, csvdelim, true));
  const auto input(loadbuf(argv[2]).second);
  if(std::strcmp(argv[1], "lword") == 0) {
    csv.insert(csv.end(), csvelim.begin(),  csvelim.end());
    csv.insert(csv.end(), csvdelim.begin(), csvdelim.end());
    const auto inputs(cutText(input, csv, delimiter));
    lword<char32_t, std::u32string> stat;
    std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> converter;
    for(int j = 0; j < inputs.size(); j ++) {
      std::u32string itrans(converter.from_bytes(inputs[j]));
      for(int i0 = 0; i0 <= input.size() / szblock; i0 ++)
        for(int i = 2; i < 20; i ++) {
          stat.init(60, i, i);
          auto words(stat.compute(itrans.substr(i0 * szblock, szblock)));
          for(auto itr = words.begin(); itr != words.end(); ++ itr)
            if(itr->str.size() > 2 && itr->count >= i) {
              std::cout << converter.to_bytes(itr->str) << ", ";
              std::cout << itr->count << std::endl;
            }
        }
    }
  } else if(std::strcmp(argv[1], "lbalance") == 0) {
    const auto cinput(cutText(input, csvelim, delimiter));
    const auto idxs(pseudoWordsBalance<double, std::string>(cinput, csv, Mbalance));
    std::cout << idxs.size() << "sets." << std::endl;
    for(int i = 0; i < idxs.size(); i ++)
      std::cout << cinput[idxs[i]] << std::endl;
  } /* else if(std::strcmp(argv[1], "corpus") == 0) {
    corpus<double, std::string> stat;
    for(int i = 0; i < input.size() / szwindow + 1; i ++) {
      stat.init(csv, 0, 120);
      const auto& words(stat.getWords());
      stat.compute(input.substr(i * szwindow, szwindow), delimiter);
      const auto& corpus(stat.getCorpus());
      std::cout << words  << std::endl;
      // std::cout << corpus << std::endl;
    }
  } */ else if(std::strcmp(argv[1], "toc") == 0 ||
            std::strcmp(argv[1], "lack") == 0) {
    std::vector<std::string> details;
    std::vector<std::string> tocs;
    std::vector<std::string> detailwords;
    std::vector<std::string> tocwords;
    bool toc(false);
    for(int iidx = 4; iidx < argc; iidx ++) {
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
    csv.insert(csv.end(), tocwords.begin(),    tocwords.end());
    csv.insert(csv.end(), detailwords.begin(), detailwords.end());
    std::sort(csv.begin(), csv.end());
    csv.erase(std::unique(csv.begin(), csv.end()), csv.end());
    std::cout << std::string("<html><head><link rel=\"stylesheet\" type=\"text/css\" href=\"../../style.css\"></head>") << std::endl;
    std::cout << std::string("<body>");
    std::cout << preparedTOC<double, std::string>(input, csv, detailwords, details, tocwords, tocs, delimiter, szwindow, 8, threshin, .125, std::strcmp(argv[1], "lack") == 0) << std::string("<hr/>") << std::endl;
    std::cout << std::string("</body></html>");
  } else if(std::strcmp(argv[1], "reconstruct") == 0) { 
    corpus<double, std::string> stat;
    stat.compute(input, delimiter, csv);
    corpushl<double, std::string> recons(stat);
    std::cout << recons.serialize();
  } else if(std::strcmp(argv[1], "redig") == 0) {
    std::vector<double> emph;
    emph.push_back(4.);
    emph.push_back(1.);
    emph.push_back(.25);
    for(int ei = 0; ei < emph.size(); ei ++) {
      for(int i = 0; i < input.size() / szwindow + 1; i ++) {
        corpus<double, std::string> stat; 
        stat.compute(input.substr(i * szwindow, szwindow), delimiter, csv);
        corpushl<double, std::string> recons(stat);
        recons.reDig(emph[ei]);
        std::cout << recons.serialize() << std::endl;
      }
      std::cout << std::endl << std::endl;
    }
  } else if(std::strcmp(argv[1], "diff") == 0) {
    std::vector<std::string> details, details2;
    std::vector<std::string> detailwords, detailwords2;
    bool second(false);
    for(int iidx = 4; iidx < argc; iidx ++) {
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
    csv.insert(csv.end(), detailwords.begin(),  detailwords.end());
    csv.insert(csv.end(), detailwords2.begin(), detailwords2.end());
    std::sort(csv.begin(), csv.end());
    csv.erase(std::unique(csv.begin(), csv.end()), csv.end());
    std::cout << std::string("<html><head><link rel=\"stylesheet\" type=\"text/css\" href=\"../../style.css\"></head>") << std::endl;
    std::cout << std::string("<body>");
    std::cout << diff<double, std::string>(input, csv, details, detailwords, details2, detailwords2, delimiter, szwindow, threshin) << std::string("<hr/>") << std::endl;
    std::cout << "</body></html>" << std::endl;
  } else if(std::strcmp(argv[1], "stat") == 0 ||
            std::strcmp(argv[1], "findroot") == 0) {
    std::vector<std::string> rdetails;
    std::vector<std::string> rdetailwords;
    for(int iidx = 4; iidx < argc; iidx ++) {
      const auto work(loadbuf(argv[iidx]));
      rdetails.push_back(work.second);
      rdetailwords.push_back(work.first);
    }
    csv.insert(csv.end(), rdetailwords.begin(), rdetailwords.end());
    std::sort(csv.begin(), csv.end());
    csv.erase(std::unique(csv.begin(), csv.end()), csv.end());
    std::cout << std::string("<html><head><link rel=\"stylesheet\" type=\"text/css\" href=\"../../style.css\"></head>") << std::endl;
    std::cout << std::string("<body>");
    std::cout << optimizeTOC<double, std::string>(input, csv, rdetails, rdetailwords, delimiter, szwindow, 8, threshin, 1., std::strcmp(argv[1], "findroot") == 0) << std::string("<hr/>") << std::endl;
    std::cout << std::string("</body></html>");
  } else if(std::strcmp(argv[1], "prep") == 0) {
    std::vector<std::string> buf;
    corpus<double, std::string> stat;
    for(int i = 0; i < input.size() / szwindow + 1; i ++) {
      stat.compute(input.substr(i * szwindow, szwindow), delimiter, csv);
      const auto work(corpushl<double,std::string>(stat).reverseLink());
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

