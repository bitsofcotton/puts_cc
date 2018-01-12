# Puts_CC
Cluster large documents from small elementary document sets and collect a TOC and make links to medium sized set of documents.
And this program is implemented with a deterministic way.

This IS NOT the program that judges it is true or false, this IS ONLY the program that seines
the score that seems to be true or that seems to be false FROM INPUT (with no deductions for logics).

Statistics definition needs to be not dense.
If it's dense, the context we get from statistics may say nothing.
So with layered dictionaries, we can calculate a little accuratery.

# Usage
    make tools
    ./tools lword < data.txt
    ./tools corpus wordlist.txt < data.txt
    ./tools toc wordlist.txt dictionaries ... -toc topics ... < data.txt
    ./tools redig wordlist.txt < data.txt
    ./tools stat wordlist.txt dictionaries ... < data.txt
    ./tools reconstruct wordlist.txt < data.txt
    ./tools diff wordlist.txt -dict ... -dict2 ... < data.txt

# Contexts
Inspired from morphological analysis softwares and parameter auto configuring algorithms' structure.   
There's preceders Tensor representation of the text set (with a way different to this) R.B. ieee 2004.

# Status
Writing optimize toc, checking and testing implementation.
And preparing for collision detecting (searching 'NOT' definition and small sized dictionaries.). 

# Demos
https://services.limpid-intensity.info/puts.php have a working sample.
Please pseudo-login with e-mail address - salt pair.

# Example of using as a library
    #include "lword.hh"
    #include "corpus.hh"
    #include <string>
    #include <iostream>
    #include <fstream>
    
    std::string input;
    // ... initialize input.
    lword<char, std::string> stat;
    stat.init(60, 2, 2, 4);
    std::vector<word_t<std::string> > words(stat.compute(input.c_str()));
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
    std::cout << cstat0.toc(tocwords, tocs, cstat0.serialize());
    
# Tips
These programs may have a algorithm that is not carefully confirmed.  
abbrev function is not ideal definition, so in fact, it is needed to be the inverse of withDetail function, but now, not so.  
serialize function is not ideal definition, so in fact, it is needed to be the inverse of corpushl(corpus(...)) but now, not so.  
This is from simply interest, so don't use this for the reason other than standing under the context.  
And this algorithm have a specification that edges remains as obscure, to get clear edge, we should use this with morphological analysis softwares. This algorithm only sees the word distances and orders that causes rough image of what is on the table.  
This program only sees what is said in the input, not the not said.  
If there's large numbers of the sets that correctly defined texts, it may optimizable out to get opposite side's claim, but this is only the what's able to be said, so this isn't what's the matter.

# Another downloads
* https://ja.osdn.net/projects/puts-cc/
* https://www.sourceforge.net/projects/puts-cc/
* https://konbu.sakura.ne.jp/files/
* https://files.limpid-intensity.info/
