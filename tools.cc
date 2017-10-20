#include <string>
#include <iostream>
#include <fstream>
#include <locale>
#include <codecvt>
#include "lword.hh"
#include "corpus.hh"
#include "corpushl.hh"
#include "file2eigen.hh"

void usage() {
  std::cout << "tools (lword|toc)" << std::endl;
}

int main(int argc, const char* argv[]) {
  std::ios::sync_with_stdio(false);
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
  while(std::getline(std::cin, buf, '\n') && !std::cin.eof() && !std::cin.bad()) input += buf + std::string("\n");
  switch(mode) {
  case 0:
    {
#if 0
      lword<char32_t, std::u32string> stat;
      std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> converter;
      std::u32string itrans(converter.from_bytes(input));
      for(int i = 2; i < 20; i ++) {
        stat.init(60, i, i);
        auto words(stat.compute(itrans.c_str()));
        for(auto itr = words.begin(); itr != words.end(); ++ itr)
          if(itr->str.size() > 2 && itr->count >= i) {
            std::cout << converter.to_bytes(itr->str) << ", ";
            std::cout << itr->count << std::endl;
          }
      }
#else
      lword<char, std::string> stat;
      for(int i = 2; i < 20; i ++) {
        stat.init(60, i, i);
        auto words(stat.compute(input.c_str()));
        for(auto itr = words.begin(); itr != words.end(); ++ itr)
          if(itr->str.size() > 2 && itr->count >= i) {
            std::cout << itr->str << ", ";
            std::cout << itr->count << std::endl;
          }
      }
#endif
      break;
    }
  case 1:
    {
      corpus<double, char> stat;
      std::string          wordsbuf, line;
      std::ifstream input2;
      input2.open(argv[2]);
      while(getline(input2, line) && !input2.eof() && !input2.bad()) wordsbuf += line + std::string("\n");
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
        std::string   line;
        std::ifstream input2;
        input2.open(argv[2]);
        while(getline(input2, line) && !input2.eof() && !input2.bad()) wordbuf += line + std::string("\n");
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
        while(getline(input2, line) && !input2.eof() && !input2.bad()) inbuf += line + std::string("\n");
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
      std::cerr << "analysing input text." << std::endl;
      for(int i = 0; i < details.size(); i ++)
        cstat0 = cstat0.withDetail(detailwords[i], details[i]);
      std::cerr << "analysing topics text." << std::endl;
      for(int i = 0; i < tocs.size(); i ++)
        for(int j = 0; j < details.size(); j ++)
          tocs[i] = tocs[i].withDetail(detailwords[j], details[j]);
      std::cerr << "getting toc." << std::endl;
      std::cout << cstat0.toc(details, detailwords, tocs, tocwords, 0, 0.);
      std::cout << std::endl << std::endl;
      corpushl<double, char> summ(cstat0);
      summ.simpleThresh(.8);
      for(int i = 0; i < tocs.size(); i ++) {
        std::vector<std::string> work(summ.reverseLink(tocs[i]));
        std::cout << tocwords[i] << " : " << std::endl;
        for(int j = 0; j < work.size(); j ++)
          std::cout << work[j] << std::endl;
      }
    }
    break;
  }
  return 0;
}

