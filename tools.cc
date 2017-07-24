#include "lword.hh"
#include "corpus.hh"
#include "corpushl.hh"
#include "file2eigen.hh"
#include <string>
#include <iostream>
#include <fstream>

void usage() {
  std::cout << "tools (lword|toc)" << std::endl;
}

int main(int argc, const char* argv[]) {
  if(argc < 2) {
    usage();
    return - 1;
  }
  int mode = - 1;
  if(std::strcmp(argv[1], "lword") == 0)
    mode = 0;
  else if(std::strcmp(argv[1], "corpus") == 0 && argc > 2)
    mode = 1;
  else if(std::strcmp(argv[1], "toc") == 0 && argc > 2)
    mode = 2;
  else {
    usage();
    return - 2;
  }
  std::string buf;
  std::string input;
  while(std::getline(std::cin, buf, '\n')) input += buf + std::string("\n");
  switch(mode) {
  case 0:
    {
      lword<char> stat;
      stat.init(120, 2, 2);
      std::vector<word_t<char> > words(stat.compute(input.c_str()));
      for(std::vector<word_t<char> >::iterator itr = words.begin(); itr != words.end(); ++ itr) {
        std::cout << itr->str << ", ";
        std::cout << itr->count << std::endl;
      }
      break;
    }
  case 1:
    {
      corpus<double, char> stat;
      std::string          wordsbuf, line;
      std::ifstream input2;
      input2.open(argv[2]);
      while(getline(input2, line)) wordsbuf += line + std::string("\n");
      input2.close();
      
      stat.init(wordsbuf.c_str(), 0, 120);
      const std::vector<std::string>& words(stat.getWords());
      stat.compute(input.c_str());
      Eigen::Matrix<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Dynamic, Eigen::Dynamic> corpus(stat.getCorpus());
      std::cout << words << std::endl;
      std::cout << corpus << std::endl;
    }
    break;
  case 2:
    {
      std::string wordbuf;
      {
        std::string   inbuf, line;
        std::ifstream input2;
        input2.open(argv[2]);
        while(getline(input2, line)) wordbuf += line + std::string("\n");
        input2.close();
      }
      std::vector<corpushl<double, char> > details;
      std::vector<corpushl<double, char> > tocs;
      std::vector<std::string>             detailwords;
      std::vector<std::string>             tocwords;
      bool toc(false);
      for(int iidx = 3; iidx < argc; iidx ++) {
        if(std::string(argv[iidx]) == std::string("-toc")) {
          toc = true;
          continue;
        }
        std::string   inbuf, line;
        std::ifstream input2;
        input2.open(argv[iidx]);
        while(getline(input2, line)) inbuf += line + std::string("\n");
        input2.close();
        corpus<double, char> cstat;
        cstat.init(wordbuf.c_str(), 0, 120);
        cstat.compute(inbuf.c_str());
        std::string wbuf(argv[iidx]);
        int slash = - 1;
        for(int j = 0; j < wbuf.size(); j ++)
          if(wbuf[j] == '/')
            slash = j;
        std::string wwbuf;
        slash ++;
        for(int j = slash; j < wbuf.size(); j ++)
          wwbuf += wbuf[j];
        if(toc) {
          tocs.push_back(corpushl<double, char>(cstat));
          tocwords.push_back(wwbuf);
        } else {
          details.push_back(corpushl<double, char>(cstat));
          detailwords.push_back(wwbuf);
        }
      }
      corpushl<double, char> cstat0;
      {
        corpus<double, char> cstat;
        cstat.init(wordbuf.c_str(), 0, 120);
        cstat.compute(input.c_str());
        cstat0 = corpushl<double, char>(cstat);
      }
      for(int i = 0; i < details.size(); i ++)
        cstat0 = cstat0.withDetail(detailwords[i], details[i]);
      std::cout << cstat0.toc(details, detailwords, tocs, tocwords, 0, 0.);
    }
    break;
  }
  return 0;
}

