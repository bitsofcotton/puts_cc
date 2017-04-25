#include "lword.hh"
#include "corpushl.hh"
#include <string>
#include <iostream>
#include <fstream>

void usage() {
  std::cout << "tools (lword|corpushl)" << std::endl;
}

int main(int argc, const char* argv[]) {
  if(argc < 1) {
    usage();
    return - 1;
  }
  int mode = - 1;
  if(std::strcmp(argv[1], "lword") == 0)
    mode = 0;
  else if(std::strcmp(argv[1], "corpushl") == 0 && argc > 2)
    mode = 1;
  else {
    usage();
    return - 2;
  }
  std::string buf;
  std::string input;
  while(std::getline(std::cin, buf, '\n')) input += buf;
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
      corpushl<double, char> stat;
      std::string   wordsbuf, line;
      std::ifstream input2;
      input2.open(argv[2]);
      while(getline(input2, line)) wordsbuf += line + '\n';
      input2.close();
      
      stat.init(wordsbuf.c_str(), 120);
      const std::vector<std::string>& words(stat.getWords());
      Eigen::Matrix<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Dynamic, Eigen::Dynamic> corpus(stat.compute(input.c_str()));
      for(int i = 0; i < corpus.rows(); i ++)
        for(int j = 0; j < corpus.cols(); j ++)
          std::cout << words[i] << " - " << words[j] << " : " << corpus(i, j).transpose() << std::endl;
    }
  }
  return 0;
}

