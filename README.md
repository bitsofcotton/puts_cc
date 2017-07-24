# Puts_CC
Cluster large documents from small elementary document sets.

This program aims to cluster and collect TOC and make links to medium sized set of documents.  
Statistics definition needs to be not dense.
If it's dense, the context we get from statistics may say nothing.
So with layered dictionaries, we can calculate a little accuratery.

# Usage
    make tools
    ./tools lword < data.txt
    ./tools corpus wordlist.txt < data.txt
    ./tools toc wordlist.txt dictionaries ... -toc topics ... < data.txt

# Example of using as a library
    #include "lword.hh"
    #include "corpus.hh"
    #include "corpushl.hh"
    #include <string>
    #include <iostream>
    #include <fstream>
    
    std::string input;
    // ... initialize input.
    lword<char> stat;
    stat.init(120, 2, 2);
    std::vector<word_t<char> > words(stat.compute(input.c_str()));
    for(auto itr = words.begin(); itr != words.end(); ++ itr) {
      std::cout << itr->str << ", ";
      std::cout << itr->count << std::endl;
    }
    
    std::string wordlist;
    // ... initialize wordlist (\t or \n separated).
    corpus<double, char> cstat;
    cstat.init(wordlist.c_str(), 0, 120);
    cstat.compute(input.c_str());
    
    corpushl<double, char> cstats;
    cstats = corpushl<double, char>(cstat);
    
    std::vector<corpushl<double, char> > details, tocs;
    std::vector<std::string> detailwords, tocwords;
    // ... initialize details and tocs.
    
    for(int i = 0; i < details.size(); i ++)
      cstat0 = cstat0.withDetail(detailwords[i], details[i]);
    // cstat0 is detailed corpus.
    
    // sample output for toc.
    std::cout << cstat0.toc(details, detailwords, details, tocwords, tocs, 10);
    
