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
  std::cout << "tools (lword|lbalance|corpus|toc|redig|stat|reconstruct|diff|getdict)" << std::endl;
}

const int szwindow(200);
const int szblock(200 * 200);
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
  if(std::strcmp(argv[1], "lword") == 0)
    mode = 0;
  else if(std::strcmp(argv[1], "lbalance") == 0 && argc > 2)
    mode = 8;
  else if(std::strcmp(argv[1], "corpus") == 0 && argc > 2)
    mode = 1;
  else if(std::strcmp(argv[1], "toc") == 0 && argc > 2)
    mode = 2;
  else if(std::strcmp(argv[1], "redig") == 0 && argc > 2)
    mode = 4;
  else if(std::strcmp(argv[1], "stat") == 0 && argc > 2)
    mode = 6;
  else if(std::strcmp(argv[1], "reconstruct") == 0 && argc > 2)
    mode = 3;
  else if(std::strcmp(argv[1], "diff") == 0 && argc > 2)
    mode = 5;
  else if(std::strcmp(argv[1], "getdict") == 0 && argc > 2)
    mode = 7;
  else {
    usage();
    return - 2;
  }
  std::string buf;
  std::string input;
  while(getline(std::cin, buf)) {
    input += buf + std::string("\n");
    if(std::cin.eof() || std::cin.bad())
      break;
  }
  switch(mode) {
  case 0:
    // lword
    {
      if(2 <= argc) {
        auto elimlist(cutText(loadbuf(argv[2]).second, csvelim, csvdelim));
        elimlist.insert(elimlist.end(), csvelim.begin(), csvelim.end());
        elimlist.push_back(std::string("\r"));
        elimlist.push_back(std::string("\n"));
        input = cutText(input, elimlist, vector<std::string>())[0];
        std::cerr << input << std::endl;
      }
      lword<char32_t, std::u32string> stat;
      std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> converter;
      std::u32string itrans(converter.from_bytes(input));
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
      break;
    }
  case 8:
    // lbalance.
    {
      std::string workbuf;
      if(2 < argc)
        workbuf = loadbuf(argv[2]).second;
      auto elims(csvelim);
      elims.insert(elims.end(), csvdelim.begin(), csvdelim.end());
      auto inputs(cutText(input, elims, delimiter));
      auto idxs(pseudoWordsBalance<double, std::string>(inputs, cutText(workbuf, csvelim, csvdelim), Mbalance));
      std::cout << idxs.size() << "sets." << std::endl;
      for(int i = 0; i < idxs.size(); i ++)
        std::cout << inputs[idxs[i]] << std::endl;
    }
    break;
  case 1:
    // corpus
    {
      corpus<double, std::string> stat;
      const auto words0(cutText(loadbuf(argv[2]).second, csvelim, csvdelim));
      for(int i = 0; i < input.size() / szwindow + 1; i ++) {
        stat.init(words0, 0, 120);
        const auto& words(stat.getWords());
        stat.compute(input.substr(i * szwindow, szwindow), delimiter);
        const auto& corpus(stat.getCorpus());
        std::cout << words  << std::endl;
        std::cout << corpus << std::endl;
      }
    }
    break;
  case 2:
    // toc
    {
      const auto words0(cutText(loadbuf(argv[2]).second, csvelim, csvdelim));
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
      std::cout << std::string("<html><head><link rel=\"stylesheet\" type=\"text/css\" href=\"../../style.css\"></head>") << std::endl;
      std::cout << std::string("<body>");
      for(int i = 0; i <= input.size() / szblock; i ++)
        std::cout << preparedTOC<double, std::string>(input.substr(i * szblock, szblock), words0, detailwords, details, tocwords, tocs, delimiter, szwindow, double(.25), .125) << std::string("<hr/>") << std::endl;
      std::cout << std::string("</body></html>");
    }
    break;
  case 3:
    // reconstruct
    {
      auto wordbuf(loadbuf(argv[2]).second);
      corpus<double, std::string> stat;
      stat.init(cutText(wordbuf, csvelim, csvdelim), 0, 120);
      stat.compute(input, delimiter);
      corpushl<double, std::string> recons(stat);
      std::cout << recons.serialize();
    }
    break;
  case 4:
    // redig
    {
      const auto words0(cutText(loadbuf(argv[2]).second, csvelim, csvdelim));
      std::vector<double> emph;
      emph.push_back(4.);
      emph.push_back(1.);
      emph.push_back(.25);
      for(int ei = 0; ei < emph.size(); ei ++) {
        for(int i = 0; i < input.size() / szwindow + 1; i ++) {
          corpus<double, std::string> stat; 
          stat.init(words0, 0, 120);
          stat.compute(input.substr(i * szwindow, szwindow), delimiter);
          corpushl<double, std::string> recons(stat);
          recons.reDig(emph[ei]);
          std::cout << recons.serialize() << std::endl;
        }
        std::cout << std::endl << std::endl;
      }
    }
    break;
  case 5:
    // diff
    {
      const auto words0(cutText(loadbuf(argv[2]).second, csvelim, csvdelim));;
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
      std::cout << std::string("<html><head><link rel=\"stylesheet\" type=\"text/css\" href=\"../../style.css\"></head>") << std::endl;
      std::cout << std::string("<body>");
      for(int i = 0; i <= input.size() / szblock; i ++)
        std::cout << diff<double, std::string>(input.substr(i * szblock, szblock), words0, details, detailwords, details2, detailwords2, delimiter, szwindow) << std::string("<hr/>") << std::endl;
      std::cout << "</body></html>" << std::endl;
    }
    break;
  case 6:
    // stat
    {
      const auto words0(cutText(loadbuf(argv[2]).second, csvelim, csvdelim));;
      std::vector<std::string> rdetails;
      std::vector<std::string> rdetailwords;
      for(int iidx = 3; iidx < argc; iidx ++) {
        const auto work(loadbuf(argv[iidx]));
        rdetails.push_back(work.second);
        rdetailwords.push_back(work.first);
      }
      std::cout << std::string("<html><head><link rel=\"stylesheet\" type=\"text/css\" href=\"../../style.css\"></head>") << std::endl;
      std::cout << std::string("<body>");
      for(int i = 0; i <= input.size() / szblock; i ++)
        std::cout << optimizeTOC<double, std::string>(input.substr(i * szblock, szblock), words0, rdetails, rdetailwords, delimiter, szwindow, 8, 1.) << std::string("<hr/>") << std::endl;
      std::cout << std::string("</body></html>");
    }
    break;
  case 7:
    // get dict.
    {
    }
    break;
  }
  return 0;
}

