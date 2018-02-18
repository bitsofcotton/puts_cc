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
  std::cout << "tools (lword|corpus|toc|redig|stat|reconstruct|diff|getdict)" << std::endl;
}

const int szwindow(200);
std::vector<std::string> delimiter;

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
//  std::ios::sync_with_stdio(false);
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
        std::string   buf;
        std::string   elimbuf;
        std::ifstream input2;
        input2.open(argv[2]);
        while(getline(input2, buf)) {
          elimbuf += buf + std::string("\n");
          if(input2.eof() || input2.bad())
            break;
        }
        input2.close();
        lword<char, std::string> eliminater;
        input = eliminater.eliminate(input, elimbuf);
        std::cerr << input << std::endl;
      }
#if 1
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
    // corpus
    {
      std::string          wordsbuf;
      std::string          line;
      std::ifstream        input2;
      input2.open(argv[2]);
      while(getline(input2, line)) {
        wordsbuf += line + std::string("\n");
        if(input2.eof() || input2.bad())
          break;
      }
      input2.close();
      corpus<double, char> stat;
      for(int i = 0; i < input.size() / szwindow + 1; i ++) {
        stat.init(wordsbuf.c_str(), 0, 120);
        const std::vector<std::string>& words(stat.getWords());
        stat.compute(input.substr(i * szwindow, szwindow).c_str(), delimiter);
        Eigen::Matrix<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Dynamic, Eigen::Dynamic> corpus(stat.getCorpus());
        std::cout << words  << std::endl;
        std::cout << corpus << std::endl;
      }
    }
    break;
  case 2:
    // toc
    {
      std::string wordbuf;
      {
        std::string   line;
        std::ifstream input2;
        input2.open(argv[2]);
        while(getline(input2, line)) {
          wordbuf += line + std::string("\n");
          if(input2.eof() || input2.bad())
            break;
        }
        input2.close();
      }
      std::vector<std::vector<corpushl<double, char> > > details;
      std::vector<std::vector<corpushl<double, char> > > tocs;
      std::vector<std::string> detailwords;
      std::vector<std::string> tocwords;
      bool toc(false);
      for(int iidx = 3; iidx < argc; iidx ++) {
        if(std::string(argv[iidx]) == std::string("-toc")) {
          toc = true;
          continue;
        }
        std::string   inbuf, line;
        std::ifstream input2;
        input2.open(argv[iidx]);
        while(getline(input2, line)) {
          inbuf += line + std::string("\n");
          if(input2.eof() || input2.bad())
            break;
        }
        input2.close();
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
          tocs.push_back(std::vector<corpushl<double, char> >());
          tocwords.push_back(wwbuf);
        } else {
          details.push_back(std::vector<corpushl<double, char> >());
          detailwords.push_back(wwbuf);
        }
        for(int i = 0; i < inbuf.size() / szwindow + 1; i ++) {
          corpus<double, char> cstat;
          cstat.init(wordbuf.c_str(), 0, 120);
          cstat.compute(inbuf.substr(i * szwindow, szwindow).c_str(), delimiter);
          if(toc)
            tocs[tocs.size() - 1].push_back(corpushl<double, char>(cstat));
          else
            details[details.size() - 1].push_back(corpushl<double, char>(cstat));
        }
      }
      std::vector<corpushl<double, char> > cstat0;
      for(int i = 0; i < input.size() / szwindow + 1; i ++) {
        corpus<double, char> cstat;
        cstat.init(wordbuf.c_str(), 0, 120);
        cstat.compute(input.substr(i * szwindow, szwindow).c_str(), delimiter);
        cstat0.push_back(corpushl<double, char>(cstat));
      }
      for(int i = 0; i < cstat0.size(); i ++)
        std::cout << cstat0[i].serialize() << std::endl;
      std::cerr << "analysing input text." << std::endl;
      for(int i = 0; i < details.size(); i ++)
        for(int j = 0; j < cstat0.size(); j ++)
          for(int k = 0; k < details[i].size(); k ++)
            cstat0[j] = cstat0[j].withDetail(detailwords[i], details[i][k]);
      std::cerr << "analysing topics text." << std::endl;
      for(int i = 0; i < tocs.size(); i ++)
        for(int j = 0; j < details.size(); j ++)
          for(int k = 0; k < tocs[i].size(); k ++)
            for(int l = 0; l < details[j].size(); l ++)
              tocs[i][k] = tocs[i][k].withDetail(detailwords[j], details[j][l]);
      std::cerr << "getting toc." << std::endl;
      for(int i = 0; i < cstat0.size(); i ++)
        std::cout << cstat0[i].toc(tocs, tocwords);
      std::cout << std::endl << std::endl;
      std::vector<corpushl<double, char> >& summ(cstat0);
      for(int ii = 0; ii < summ.size(); ii ++) {
        summ[ii] = summ[ii].simpleThresh(.8);
        for(int i = 0; i < tocs.size(); i ++) {
          std::cout << ii << " - " <<  tocwords[i] << " : " << std::endl;
          std::vector<std::string> workb;
          for(int j = 0; j < tocs[i].size(); j ++) {
            std::vector<std::string> work(summ[ii].reverseLink(tocs[i][j]));
            if(work != workb) {
              for(int k = 0; k < work.size(); k ++)
                std::cout << work[k] << ", ";
              std::cout << std::endl;
              workb = work;
            }
          }
        }
      }
    }
    break;
  case 3:
    // reconstruct
    {
      std::string wordbuf;
      {
        std::string   line;
        std::ifstream input2;
        input2.open(argv[2]);
        while(getline(input2, line)) {
          wordbuf += line + std::string("\n");
          if(input2.eof() || input2.bad())
            break;
        }
        input2.close();
      }
      corpus<double, char> stat;
      stat.init(wordbuf.c_str(), 0, 120);
      stat.compute(input.c_str(), delimiter);
      corpushl<double, char> recons(stat);
      std::cout << recons.serialize();
    }
    break;
  case 4:
    // redig
    {
      std::string wordbuf;
      { 
        std::string   line;
        std::ifstream input2;
        input2.open(argv[2]); 
        while(getline(input2, line)) {
          wordbuf += line + std::string("\n");
          if(input2.eof() || input2.bad())
            break;
        }
        input2.close();
      }
      std::vector<double> emph;
      emph.push_back(4.);
      emph.push_back(1.);
      emph.push_back(.25);
      for(int ei = 0; ei < emph.size(); ei ++) {
        for(int i = 0; i < input.size() / szwindow + 1; i ++) {
          corpus<double, char> stat; 
          stat.init(wordbuf.c_str(), 0, 120);
          stat.compute(input.substr(i * szwindow, szwindow).c_str(), delimiter);
          corpushl<double, char> recons(stat);
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
      std::string wordbuf;
      {
        std::string   line;
        std::ifstream input2;
        input2.open(argv[2]);
        while(getline(input2, line)) {
          wordbuf += line + std::string("\n");
          if(input2.eof() || input2.bad())
            break;
        }
        input2.close();
      }
      std::vector<corpushl<double, char> > details, details2;
      std::vector<std::string>             detailwords, detailwords2;
      bool second(false);
      for(int iidx = 3; iidx < argc; iidx ++) {
        if(std::string(argv[iidx]) == std::string("-dict")) {
          second = false;
          continue;
        } else if(std::string(argv[iidx]) == std::string("-dict2")) {
          second = true;
          continue;
        }
        std::string   inbuf, line;
        std::ifstream input2;
        input2.open(argv[iidx]);
        while(getline(input2, line)) {
          inbuf += line + std::string("\n");
          if(input2.eof() || input2.bad())
            break;
        }
        input2.close();
        corpus<double, char> cstat;
        cstat.init(wordbuf.c_str(), 0, 120);
        cstat.compute(inbuf.c_str(), delimiter);
        std::string wbuf(argv[iidx]);
        int slash(- 1);
        for(int j = 0; j < wbuf.size(); j ++)
          if(wbuf[j] == '/')
            slash = j;
        std::string wwbuf;
        slash ++;
        for(int j = slash; j < wbuf.size(); j ++)
          wwbuf += wbuf[j];
        if(second) {
          details2.push_back(corpushl<double, char>(cstat));
          detailwords2.push_back(wwbuf);
        } else {
          details.push_back(corpushl<double, char>(cstat));
          detailwords.push_back(wwbuf);
        }
      }
      std::cerr << "analysing input text." << std::endl;
      std::cout << std::string("<html><head><link rel=\"stylesheet\" type=\"text/css\" href=\"../../style.css\"></head>") << std::endl;
      std::cout << std::string("<body>");
      for(int i = 0; i < input.size() / szwindow + 1; i ++) {
        corpus<double, char> cstatorig;
        cstatorig.init(wordbuf.c_str(), 0, 120);
        cstatorig.compute(input.substr(i * szwindow, szwindow).c_str(), delimiter);
        corpushl<double, char> cstat0(cstatorig), cstat1(cstatorig);
        for(int i = 0; i < details.size(); i ++)
          cstat0 = cstat0.withDetail(detailwords[i] , details[i]);
        cstat0.reDig(double(4));
        for(int i = 0; i < details2.size(); i ++)
          cstat1 = cstat1.withDetail(detailwords2[i], details2[i]);
        cstat1.reDig(double(4));
        if(cstat0 != cstat1) {
          auto diff(cstat0 - cstat1);
          diff = diff.simpleThresh(1);
          std::cout << diff.serialize() << "<br/>" << std::endl;
          std::cout << diff.reverseLink(cstatorig) << std::endl;
          std::cout << "<br/><br/>" << std::endl;
        }
      }
      std::cout << "</body></html>" << std::endl;
    }
    break;
  case 6:
    // stat
    {
      std::string wordbuf;
      {
        std::string   line;
        std::ifstream input2;
        input2.open(argv[2]);
        while(getline(input2, line)) {
          wordbuf += line + std::string("\n");
          if(input2.eof() || input2.bad())
            break;
        }
        input2.close();
      }
      std::vector<std::string> rdetails;
      std::vector<std::string> rdetailwords;
      for(int iidx = 3; iidx < argc; iidx ++) {
        std::string   inbuf, line;
        std::ifstream input2;
        input2.open(argv[iidx]);
        while(getline(input2, line)) {
          inbuf += line + std::string("\n");
          if(input2.eof() || input2.bad())
            break;
        }
        input2.close();
        std::string wbuf(argv[iidx]);
        int slash = - 1;
        for(int j = 0; j < wbuf.size(); j ++)
          if(wbuf[j] == '/')
            slash = j;
        std::string wwbuf;
        slash ++;
        for(int j = slash; j < wbuf.size(); j ++)
          wwbuf += wbuf[j];
        rdetails.push_back(inbuf);
        rdetailwords.push_back(wwbuf);
      }
      const int tot_cont(max(min(12, int(input.size()) / szwindow), 1));
      std::cout << std::string("<html><head><link rel=\"stylesheet\" type=\"text/css\" href=\"../../style.css\"></head>") << std::endl;
      std::cout << std::string("<body>");
      std::cout << optimizeTOC<double, char>(input, wordbuf.c_str(), rdetails, rdetailwords, delimiter, szwindow, tot_cont, 8, .125) << std::endl;
      std::cout << optimizeTOC<double, char>(input, wordbuf.c_str(), rdetails, rdetailwords, delimiter, szwindow, tot_cont, 8, 1.) << std::endl;
      std::cout << optimizeTOC<double, char>(input, wordbuf.c_str(), rdetails, rdetailwords, delimiter, szwindow, tot_cont, 8, 8.) << std::endl;
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

