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
    ./tools corpushl wordlist.txt < data.txt
