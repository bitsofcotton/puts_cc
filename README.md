# Puts_CC
Cluster large documents from small elementary document sets.

This program aims to cluster and collect TOC and make links to medium sized set of documents.
So first, collect unknown words, then, make a word distance table using prepared word list and collected ones.
Second, make a detail on that with using prepared word distance tables, then,
make a TOC and links from statistical information.  
And using detailed statistics, we may able to detect collisions or which prejudice gains with prejudice rule set.

Now, only collecting words and making word distance table is available.

# Usage
    make tools
    ./tools lword < data.txt
    ./tools corpus wordlist.txt < data.txt
    ./tools detail wordlist.txt dictionaries ... < data.txt

# Example of using as a library
    #include "lword.hh"
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
    
    std::vector<corpushl<double, char> > details;
    std::vector<std::string> detailwords;
    // initialize details.
    for(int i = 0; i < details.size(); i ++)
      cstat0 = cstat0.withDetail(detailwords[i], details[i]);
    // cstat0 is detailed corpus.
    
    // sample output for toc.
    std::cout << cstat0.toc(details, detailwords, details, 10);
    
