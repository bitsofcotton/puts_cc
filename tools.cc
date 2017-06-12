#include "lword.hh"
#include "corpus.hh"
#include "corpushl.hh"
#include "file2eigen.hh"
#include <string>
#include <iostream>
#include <fstream>

void usage() {
  std::cout << "tools (lword|corpus|corpusp|lcp)" << std::endl;
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
  else if(std::strcmp(argv[1], "corpusp") == 0 && argc > 3)
    mode = 2;
  else if(std::strcmp(argv[1], "lcp") == 0 && argc > 2)
    mode = 3;
  else if(std::strcmp(argv[1], "detail") == 0 && argc > 2)
    mode = 4;
  else {
    usage();
    return - 2;
  }
  std::string buf;
  std::string input;
  while(std::getline(std::cin, buf, '\n')) input += buf + "\n";
  switch(mode) {
  case 0:
    {
      lword<char> stat;
      stat.init(120, 2, 2);
      std::vector<word_t<char> > words(stat.compute(input.c_str()));
      for(auto itr = words.begin(); itr != words.end(); ++ itr) {
        std::cout << itr->str << ", ";
        std::cout << itr->count << std::endl;
      }
      break;
    }
  case 1:
    {
      corpus<double, char> stat;
      std::string   wordsbuf, line;
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
      corpus<double, char> stat, stat2;
      std::string wordsbuf, line;
      std::ifstream input2;
      input2.open(argv[2]);
      while(getline(input2, line)) wordsbuf += line + std::string("\n");
      input2.close();
      
      std::string wordsbuf2;
      input2.open(argv[3]);
      while(getline(input2, line)) wordsbuf2 += line + std::string("\n");
      input2.close();
      
      stat.init(wordsbuf.c_str(), 0, 120);
      const std::vector<std::string>& words(stat.getWords());
      stat.compute(input.c_str());
      stat2.init(wordsbuf.c_str(), 0, 120);
      const std::vector<std::string>& words2(stat2.getWords());
      stat2.compute(wordsbuf2.c_str());
      corpushl<double, char> hlstat(stat), hlstat2(stat2), res;
      res = hlstat + hlstat2;
      Eigen::Matrix<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Dynamic, Eigen::Dynamic> corpus(res.getCorpus());
      std::cout << words << std::endl;
      std::cout << corpus << std::endl;
    }
    break;
  case 3:
    {
      corpushl<double, char> last;
      for(int iidx = 2; iidx < argc; iidx ++) {
        std::string   inbuf, line;
        std::ifstream input2;
        input2.open(argv[iidx]);
        while(getline(input2, line)) inbuf += line + std::string("\n");
        input2.close();
        lword<char> stat;
        stat.init(120, 2, 2);
        std::vector<word_t<char> > words(stat.compute(inbuf.c_str()));
        for(auto itr = words.begin(); itr != words.end(); ++ itr)
          if(itr->count > 2)
            input += itr->str + std::string("\n");
        corpus<double, char> cstat;
        cstat.init(input.c_str(), 0, 120);
        cstat.compute(inbuf.c_str());
        last += corpushl<double, char>(cstat);
      }
      std::cout << last.getWords() << std::endl;
      // std::cout << last.getCorpus() << std::endl;
    }
    break;
  case 4:
    {
      std::vector<corpushl<double, char> > details;
      std::vector<std::string>             detailwords;
      for(int iidx = 3; iidx < argc; iidx ++) {
        std::string   inbuf, line;
        std::ifstream input2;
        input2.open(argv[iidx]);
        while(getline(input2, line)) inbuf += line + std::string("\n");
        input2.close();
        corpus<double, char> cstat;
        cstat.init(input.c_str(), 0, 120);
        cstat.compute(inbuf.c_str());
        details.push_back(corpushl<double, char>(cstat));
        detailwords.push_back(std::string(argv[iidx]));
      }
      corpushl<double, char> cstat0;
      {
        std::string   inbuf, line;
        std::ifstream input2;
        input2.open(argv[2]);
        while(getline(input2, line)) inbuf += line + std::string("\n");
        input2.close();
        corpus<double, char> cstat;
        cstat.init(input.c_str(), 0, 120);
        cstat.compute(inbuf.c_str());
        std::cout << cstat.getWords() << std::endl;
        cstat0 = corpushl<double, char>(cstat);
      }
      for(int i = 0; i < details.size(); i ++)
        cstat0 = cstat0.withDetail(detailwords[i], details[i]);
      std::cout << cstat0.getWords() << std::endl;
      // std::cout << cstat0.getCorpus() << std::endl;
    }
    break;
  }
  return 0;
}

