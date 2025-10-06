# Puts_CC
Cluster large documents from small elementary document sets and collect a TOC and make links to medium sized set of documents.

This IS NOT the program judges it is true or false, this IS ONLY the program seins the score it seems to be true or it seems to be false FROM INPUT (with NO deductions for logics).

# Usage
    make puts
    ./puts lword prepared_word_list.txt? < data.txt
    ./puts lbalance wordlist.txt? < data.txt
    ./puts toc   wordlist.txt dictionaries ... -toc topics ... < data.txt
    ./puts lack  wordlist.txt dictionaries ... -toc topics ... < data.txt
    ./puts redig wordlist.txt? < data.txt
    ./puts stat (wordlist.txt? (dictionaries ...)?)? < data.txt
    ./puts findroot (wordlist.txt? (dictionaries ...)?)? < data.txt
    ./puts diff  wordlist.txt? -dict ... -dict2 ... < data.txt
    ./puts same  wordlist.txt? -dict ... -dict2 ... < data.txt
    ./puts prep  wordlist.txt? < data.txt
    ./puts pred (wordlist.txt? (dictionaries ...)?)? < data.txt

# Contexts
Inspired from morphological analysis softwares and parameter auto configuring algorithms' structure.   
There's preceders Tensor representation of the text set (with a way different to this) R.B. ieee 2004.

# Tips
Implanted tips into lieonn.hh as a comment. Also analyse.php as a comment.

# Another downloads
* https://konbu.azurewebsites.net/ (Sample Site)
* https://drive.google.com/drive/folders/1B71X1BMttL6yyi76REeOTNRrpopO8EAR?usp=sharing
* https://1drv.ms/u/s!AnqkwcwMjB_PaDIfXya_M3-aLXw?e=qzfKcU
* https://ja.osdn.net/users/bitsofcotton/

# Archived
This repository is archived, so without bug report, will no change.
2023/03/13 integrate some files into lieonn.hh.
2023/03/24 code clean.
2023/07/15 add fix predTOC.
2023/09/12 add and fix predvResizeSTens.
2023/09/19 merge latest lieonn causes predv correction.
2023/09/25 change prediction strategy, prefer to use complement.
2023/10/22 we select the speed instead of accuracy on prediction, so we should to shrink the output tensor.
2023/10/29 pred output fix, also have serializeSub function additional large effect fix. Either we have prediction noise glitch which will be reduced sqrt scale tensor but we don't implement them because of huge memory usage.
2023/10/30 copy structure reliably with randtools meaning.
2023/12/15 merge lieonn.hh causes pred takes invariant only once in the whole.
2023/12/17 serializeSub minus tensor to be reverse order of abs(minus tensor).
2023/01/20 Auto configure szwindow with sqrt input size.
2023/01/21 Reconfigure default scorethresh.
2024/04/02 P01 fix.
2024/04/04 only use large accuracy on calculating pnextcache, but this is broken with cache naming.
2024/04/09 merge latest lieonn from ddpmopt.
2024/04/14 merge latest lieonn, no logic change.
2024/05/17 omit prediction output.
2024/05/19 add prednoword command.
2024/05/20 reconfigure param nrwords, it's still large.
2024/05/21 fix opt direction. they doesn't need such large mem. code cleaning. improved predSTen mem usage.
2024/06/07 merge latest lieonn.
2024/06/18 merge latest ddpmopt results.
2024/06/19 merge latest ddpmopt results.
2024/06/21 merge latest ddpmopt results.
2024/06/22 p01 fatal fix. make/revert program invariant change friendly to predictors.
2024/06/23 merge latest ddpmopt.
2024/06/23 merge latest ddpmopt.
2024/06/24 merge latest ddpmopt.
2024/06/26 merge latest ddpmopt.
2024/06/26 merge latest ddpmopt.
2024/06/27 merge latest ddpmopt. predSTen bugs might be fixed.
2024/06/30 merge latest ddpmopt, change pred output to raw text with predicted and original.
2024/07/06 merge latest ddpmopt.
2024/07/07 merge latest ddpmopt.
2024/07/08 merge latest ddpmopt.
2024/07/09 predSTen code clean.
2024/08/08 update readme.
2024/08/09 fix typo.
2024/08/18 merge update from ddpmopt.
2024/09/28 merge latest ddpmopt result.
2024/10/26 predSTen output fix on counting iteration causes only 1-indexed output.
2024/12/06 merge latest ddpmopt result causes prediction improve.
2024/12/26 merge latest ddpmopt result causes prediction improve.
2025/01/27 merge latest lieonn causes improve memory usage need half ones differed to before.
2025/02/05 merge latest lieonn affects pred cmd.
2025/03/05 merge latest lieonn, avoiding one of the dimension being attacked condition.
2025/03/09 merge latest lieonn.
2025/03/17 normalize predTOC before to serialize.
2025/04/17 merge latest dimension auto tuner from ddpmopt.
2025/04/18 step parameter strategy change, no logic change on this binary.
2025/04/19 merge latest lieonn.
2025/05/16 merge latest lieonn from ddpmopt result.
2025/05/23 merge latest ddpmopt result as feed_much.
2025/05/25 merge latest ddpmopt result.
2025/06/08 merge latest ddpmopt result.
2025/06/12 compat compile option with one variant of gcc2.95.3.
2025/06/17 merge latest ddpmopt fix.
2025/06/19 merge latest ddpmopt fix.
2025/06/22 merge latest ddpmopt result.
2025/06/23 code clean, flush.
2025/06/25 merge latest ddpmopt result. also need to implant readme.md into lieonn.hh.
2025/06/25 readme.md cleaning and implant into lieonn.hh, analyse.php comment.
2025/06/28 refactor and fix around lieonn, re-compat with gcc2953.
2025/06/29-30 merge latest ddpmopt result causes pred change.
2025/07/01 merge latest ddpmopt result causes various bug fixes.
2025/07/02-03 merge latest ddpmopt result causes new prediction algorithm debug ok.
2025/07/04 merge latest ddpmopt debug ok.
2025/07/06 brush up lieonn. bug fixes. revertByProgramInvariant important fix.
2025/07/10 merge latest ddpmopt result.
2025/07/13 merge latest ddpmopt result.
2025/07/14-16 merge latest ddpmopt result, however, predSTen argments isn't optimal we have now.
2025/07/17-19 merge latest ddpmopt result, either predSTen argments isn't optimal we have now.
2025/07/20 merge latest ddpmopt result. configured predSTen params.
2025/07/24 merge latest ddpmopt result. hard configured predSTen params.
2025/07/25 merge latest ddpmopt result.
2025/07/26-28 merge latest ddpmopt result.
2025/07/28 merge latest ddpmopt result.
2025/08/01 merge latest ddpmopt result.
2025/08/02-03 merge latest ddpmopt result.
2025/08/04-06 merge latest p2 result.
2025/08/11 merge latest p2 result.
2025/08/12-15 merge latest p2 result.
2025/08/16 merge latest p2 result.
2025/08/16 merge latest p2 result.
2025/08/17-23 merge latest p2 result.
2025/08/25 merge latest ddpmopt result also fix crash.
2025/09/01 offset bug fix on pAppendMeasure.
2025/09/05 fix to compile with latest g++, pAppendMeasure fix (should be done).
2025/09/25 merge latest ddpmopt result cause bug fix around pRS, p01next.
2025/10/06 merge latest ddpmopt result cause intended fix.

