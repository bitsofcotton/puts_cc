# Puts_CC
Cluster large documents from small elementary document sets and collect a TOC and make links to medium sized set of documents.
And this program is implemented with a deterministic way. So don't use machine learning methods.

This IS NOT the program that judges it is true or false, this IS ONLY the program that seines
the score that seems to be true or that seems to be false FROM INPUT (with NO deductions for logics).

This program uses large amount of region on the memory. And, very slow because of the scans on large regions. And, almost all of memory caches are not works on this.

# Usage
    make tools
    ./tools lword < data.txt
    ./tools lword prepared_word_list.txt < data.txt
    ./tools lbalance wordlist.txt < data.txt
    ./tools toc  wordlist.txt dictionaries ... -toc topics ... < data.txt
    ./tools lack wordlist.txt dictionaries ... -toc topics ... < data.txt
    ./tools redig wordlist.txt < data.txt
    ./tools stat wordlist.txt dictionaries ... < data.txt
    ./tools findroot wordlist.txt dictionaries ... < data.txt
    ./tools diff wordlist.txt -dict ... -dict2 ... < data.txt
    ./tools same wordlist.txt -dict ... -dict2 ... < data.txt
    ./tools prep wordlist.txt < data.txt
    ./tools pred wordlist.txt (-|db.txt) dict1 ... < data.txt

# Memory Fragments
There exists better efficient memory management programs like jemalloc.
So with them, this program memory heap usage dramatically increases.

# Contexts
Inspired from morphological analysis softwares and parameter auto configuring algorithms' structure.   
There's preceders Tensor representation of the text set (with a way different to this) R.B. ieee 2004.

# Status
Checking the implementation before to freeze.
abbrev function is a little unstable now.

# Example of using as a library
Please refer tools.cc, and please include with namespace block but include guard should harms.

# Example of using as web interface
Please move the files to directory that php files and prep.py file and tools executable as puts executable.
And please mkdir data directory as configured permissions.  
And for normal use, mysql compatible database and php executable environment and python-mysql environment is needed.
And, please configure the executable paths.

prep.py needs legacy \_mysql library with python2.7, this can cause compatibility problem, so then, please rewrite prep.py file. And, database references are done only with prep.py, please configure first.

# Importing wiktionary database (GFDL)
We needs mediawiki and xml2sql softwares. And if you need, please import wiktionary huge database into your database to run with. And, dictating wiktionary needs huge amounts of memory.

# Tips
These programs may have a algorithm that is *not carefully confirmed*.  
Statistics definition needs to be not dense. If it's dense, the context we get from statistics may say nothing. So with layered dictionaries, we can calculate a little accuratery.  
And this algorithm have a specification that edges remains as obscure, to get clear edges, we should use this with morphological analysis softwares. This algorithm only sees the word distances and orders that causes rough image of what is on the table.  
This program only sees what is said in the input, not the not said.  
abbrev is needed to be the inverse of withDetail function, but now, NOT so.  
serialize function is needed to be the inverse of corpushl(corpus(...)) but now, NOT so.  
This is from simply interest, so don't use this for the reason other than standing under the context.  
If there's large numbers of the sets that correctly defined texts, it may optimizable out to get opposite side's claim, but this is only the what's able to be said, so this isn't what's the matter.  
If the dictionary we use is not valid for the text input, we cannot make text orders correct, this is because we calculate text orders by the dictionaries only.  

# Another downloads
* https://konbu.azurewebsites.net/ (Sample Site)
* https://drive.google.com/drive/folders/1B71X1BMttL6yyi76REeOTNRrpopO8EAR?usp=sharing
* https://1drv.ms/u/s!AnqkwcwMjB_PaDIfXya_M3-aLXw?e=qzfKcU
* https://ja.osdn.net/users/bitsofcotton/

# The things undone
If we use similar things with category - relation - region of measure tensor, and with morphological analysis, we might get better results and low memory condition, but there's NO implementation here. So that isn't needed by making TOCs. And if so, if there's two or more phenomenon or people, shrink into one table, it's meaningless. In such case, this program returns broken results.  
And, if we use this program to produce brand-new data from input, we need some materials and/or relations to view from and on, and so on. And, if we use this for produce some kind of stories, we need multiple tables (on the destination and on the source of identity) transition, and some relation on the world that to be described, so to apply this for the kind of it, entity space and glues is needed other than dictating.  
And, if we try to use this program to produce some mob character want to say, we need what tend to been saw by the character, this can be made a role by the sets of dictionary but, if the characters' state is unchanged, we can't make conversation with.

If we use NOT word with accurate handle, we can roughly (but deeply differed to now) divide some insists.

# Archived
This repository is archived, so without bug report, will no change.
2023/03/13 integrate some files into lieonn.hh.
2023/03/24 code clean.
2023/07/15 add fix predTOC.
2023/09/12 add and fix predvResizeSTens.
2023/09/19 merge latest lieonn causes predv correction.
2023/09/25 change prediction strategy, prefer to use complement.

