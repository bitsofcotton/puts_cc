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
  std::cout << "tools (lword|corpus|toc|redig|stat|reconstruct|diff)" << std::endl;
}

const int szwindow(300);
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
    {
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
        std::cout << cstat0[i].serialize(.9, .01) << std::endl;
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
      for(int j = 0; j < tocs.size(); j ++)
        for(int i = 0; i < cstat0.size(); i ++)
          std::cout << cstat0[i].toc(details, detailwords, tocs[j], tocwords[j], 0, .8);
      std::cout << std::endl << std::endl;
      std::vector<corpushl<double, char> > summ(cstat0);
      for(int i = 0; i < summ.size(); i ++)
        summ[i].simpleThresh(.8);
      for(int i = 0; i < tocs.size(); i ++)
        for(int j = 0; j < tocs[i].size(); j ++) {
          std::vector<std::string> work(summ[i].reverseLink(tocs[i][j]));
          std::cout << tocwords[i] << " : " << std::endl;
          for(int k = 0; k < work.size(); k ++)
            std::cout << work[k] << std::endl;
        }
    }
    break;
  case 3:
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
      std::cout << recons.serialize(.9, .01);
    }
    break;
  case 4:
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
          std::cout << recons.serialize(.9, pow(.01, emph[ei])) << std::endl;
        }
        std::cout << std::endl << std::endl;
      }
    }
    break;
  case 5:
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
        int slash = - 1;
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
      corpushl<double, char> cstat0;
      {
        corpus<double, char> cstat;
        cstat.init(wordbuf.c_str(), 0, 120);
        cstat.compute(input.c_str(), delimiter);
        cstat0 = corpushl<double, char>(cstat);
      }
      corpushl<double, char> cstat1(cstat0), cstat2(cstat0);
      std::cout << cstat0.serialize(.9, .01, .001) << std::endl;
      std::cerr << "analysing input text." << std::endl;
      for(int i = 0; i < details.size(); i ++)
        cstat0 = cstat0.withDetail(detailwords[i] , details[i]);
      for(int i = 0; i < details2.size(); i ++)
        cstat1 = cstat1.withDetail(detailwords2[i], details2[i]);
      cstat0.reDig(double(4));
      cstat1.reDig(double(4));
      auto diff(cstat0 - cstat1);
      std::cout << diff.serialize(.9, .01, .001) << std::endl;
      std::cout << cstat2.toc(details2, detailwords2, std::vector<corpushl<double, char> >(), std::string(), 0, .8) << std::endl;
    }
    break;
  case 6:
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
      std::vector<std::string> detailwords;
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
        details.push_back(std::vector<corpushl<double, char> >());
        detailwords.push_back(wwbuf);
        for(int i = 0; i < inbuf.size() / szwindow + 1; i ++) {
          corpus<double, char> cstat;
          cstat.init(wordbuf.c_str(), 0, 120);
          cstat.compute(inbuf.substr(i * szwindow, szwindow).c_str(), delimiter);
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
      std::cerr << "analysing input text." << std::endl;
      for(int i = 0; i < details.size(); i ++)
        for(int j = 0; j < cstat0.size(); j ++)
          for(int k = 0; k < details[i].size(); k ++)
            cstat0[j] = cstat0[j].withDetail(detailwords[i], details[i][k]);
      std::vector<std::vector<std::pair<double, int> > > cstats;
      for(int i = 0; i < cstat0.size(); i ++) {
        cstats.push_back(std::vector<std::pair<double, int> >());
        for(int j = 0; j < i; j ++)
          cstats[i].push_back(cstats[j][i]);
        for(int j = i; j < cstat0.size(); j ++) {
          std::cerr << "." << std::flush;
          cstats[i].push_back(std::make_pair(cstat0[i].distanceInUnion(cstat0[j]), j));
        }
      }
      for(int i = 0; i < cstats.size(); i ++)
        std::sort(cstats[i].begin(), cstats[i].end());
      std::vector<int> phrases;
      const int depth(6);
      for(int i = 0; i < depth; i ++) {
        std::vector<std::pair<double, int> > work;
        std::vector<std::vector<int> >       idxs;
        for(int j = 0; j < cstats.size(); j ++) {
          double score(0);
          idxs.push_back(std::vector<int>());
          for(int k = 0, kk = 0; k < cstats.size() && kk < cstats.size() / depth; k ++) {
            const int k0(cstats.size() - k - 1);
            if(!std::binary_search(phrases.begin(), phrases.end(), cstats[j][k0].second)) {
              score += cstats[j][k0].first;
              idxs[j].push_back(cstats[j][k0].second);
              kk ++;
            }
          }
          work.push_back(std::make_pair(score, j));
        }
        sort(work.begin(), work.end());
        const int jj(work[work.size() - 1].second);
        corpushl<double, char> cs(cstat0[jj]);
        for(int j = 0; j < idxs[jj].size(); j ++) {
          phrases.push_back(idxs[jj][idxs[jj].size() - j - 1]);
          cs += cstat0[idxs[jj][idxs[jj].size() - j - 1]];
        }
        cs *= (double(1) / idxs[jj].size());
        std::cout << cs.serialize(.75, .00001) << std::endl;
        std::sort(phrases.begin(), phrases.end());
        phrases.erase(std::unique(phrases.begin(), phrases.end()), phrases.end());
      }
    }
    break;
  }
  return 0;
}

