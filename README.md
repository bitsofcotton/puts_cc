# Puts_CC
Cluster large documents from small elementary document sets and collect a TOC and make links to medium sized set of documents.
And this program is implemented with a deterministic way.

This IS NOT the program that judges it is true or false, this IS ONLY the program that seines
the score that seems to be true or that seems to be false FROM INPUT (with NO deductions for logics).

This program uses large amount of region on the memory. And, very slow because of the scans on large regions. And, almost all of memory caches are not works on this.

# Usage
    make tools
    ./tools lword < data.txt
    ./tools lword prepared_word_list.txt < data.txt
    ./tools lbalance wordlist.txt < data.txt
    ./tools corpus wordlist.txt < data.txt
    ./tools toc wordlist.txt dictionaries ... -toc topics ... < data.txt
    ./tools redig wordlist.txt < data.txt
    ./tools stat wordlist.txt dictionaries ... < data.txt
    ./tools reconstruct wordlist.txt < data.txt
    ./tools diff wordlist.txt -dict ... -dict2 ... < data.txt
    ./tools prep wordlist.txt < data.txt

# Memory Fragments
There exists better efficient memory management programs like jemalloc.
So with them, this program memory heap usage dramatically increases.

# Contexts
Inspired from morphological analysis softwares and parameter auto configuring algorithms' structure.   
There's preceders Tensor representation of the text set (with a way different to this) R.B. ieee 2004.

# Status
Checking and testing implementation through test.cc . lword.hh, corpus.hh is a little stable, others are being in test.

# Demos
https://services.limpid-intensity.info/puts.php have a working sample.
Please pseudo-login with e-mail address - salt pair.

# Example of using as a library
Please refer tools.cc, and please include with namespace block but include guard should harms.

# Example of using as web interface
Please move to directory that php file and prep.py file and tools executable as puts executable.
And please mkdir data directory as configured permissions.  
And for normal use, mysql compatible database and php executable environment and python-mysql environment is needed.
And, please configure the executable paths.

prep.py needs legacy \_mysql library with python2.7, this can cause compatibility problem, so then, please rewrite prep.py file. And, database references are done only with prep.py, please configure first.

# Importing wiktionary database (GFDL)
We needs mediawiki and xml2sql softwares. And if you need, please import wiktionary huge database into your database to run with.

# Tips
These programs may have a algorithm that is *not carefully confirmed*.  
Statistics definition needs to be not dense. If it's dense, the context we get from statistics may say nothing. So with layered dictionaries, we can calculate a little accuratery.  
And this algorithm have a specification that edges remains as obscure, to get clear edges, we should use this with morphological analysis softwares. This algorithm only sees the word distances and orders that causes rough image of what is on the table.  
This program only sees what is said in the input, not the not said.  
abbrev is needed to be the inverse of withDetail function, but now, NOT so.  
serialize function is needed to be the inverse of corpushl(corpus(...)) but now, NOT so.  
This is from simply interest, so don't use this for the reason other than standing under the context.  
If there's large numbers of the sets that correctly defined texts, it may optimizable out to get opposite side's claim, but this is only the what's able to be said, so this isn't what's the matter.  

# Another downloads
* https://ja.osdn.net/projects/puts-cc/
* https://www.sourceforge.net/projects/puts-cc/
* https://files.limpid-intensity.info/

# The things undone
If we use similar things with category - relation - region of measure tensor, and with morphological analysis, we might get better results and low memory condition, but there's NO implementation here. So that isn't needed by making TOCs. And if so, if there's two or more phenomenon or people, shrink into one table, it's meaningless. In such case, this program returns broken results.

# The things didn't determined.
If we use NOT word with accurate handle, we can roughly (but deeply differed to now) divide some insists.
Being determined whether implement or not.
