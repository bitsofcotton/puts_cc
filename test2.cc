#include <string>
#include <iterator>
#include <iostream>
#include <fstream>
#include <locale>
#include <codecvt>
#include "lword.hh"
#include "corpus.hh"
#include "file2eigen.hh"

void usage() {
  std::cout << "tools stat" << std::endl;
}

const int szwindow(200);
const int szblock(20000);
const int Mbalance(40);
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
  int mode = - 1;
  if(argc < 2) {
    usage();
    return - 1;
  }
  if(std::strcmp(argv[1], "stat") == 0 && argc > 2)
    mode = 6;
  else {
    usage();
    return - 2;
  }
  switch(mode) {
  case 6:
    // stat : statistics with optimized TOC with link to original articles.
    {
      const auto words0(cutText(loadbuf(argv[2]).second, csvelim, csvdelim));;
      const auto input(loadbuf(argv[3]).second);
      std::vector<std::string> rdetails;
      std::vector<std::string> rdetailwords;
      for(int iidx = 4; iidx < argc; iidx ++) {
        const auto work(loadbuf(argv[iidx]));
        rdetails.push_back(work.second);
        rdetailwords.push_back(work.first);
      }
      std::cout << std::string("<html><head><link rel=\"stylesheet\" type=\"text/css\" href=\"../../style.css\"></head>") << std::endl;
      std::cout << std::string("<body>");
      for(int i = 0; i <= input.size() / szblock; i ++)
        std::cout << optimizeTOC<double, std::string>(input.substr(i * szblock, szblock), std::string("ref") + std::to_string(i) + std::string("-"), words0, rdetails, rdetailwords, delimiter, szwindow, 8, 1.) << std::string("<hr/>") << std::endl;
      std::cout << std::string("</body></html>");
    }
    break;
  }
  return 0;
}

