# Puts_CC
中程度のテキストを初等的な辞書を使って個々の集合に分けるプログラムです。

このプログラムは、個々の集合に分けた後に、目次と元記事へのリンクを作成することを目的としています。  
また、このプログラムは入力された情報に対して何ら真偽を判定するものではなく、入力された情報から、
それが真と主張されているか、または偽と主張されているか、を表記からスコア付けするものです。  
内部処理で使用されるテンソルが Dense な場合、示される情報は意味を失っています。
複数レイヤに別れた辞書をうまく使うと、このような状態を避けることができます。

# 使い方
    make tools
    ./tools lword < data.txt
    ./tools corpus wordlist.txt < data.txt
    ./tools toc wordlist.txt dictionaries ... -toc topics ... < data.txt
    ./tools redig wordlist.txt < data.txt
    ./tools reconstruct wordlist.txt < data.txt
    ./tools diff wordlist.txt -dict ... -dict2 ... < data.txt

# 文脈
形態要素解析のソフトウェア、及び自動パラメタ調整のソフトウェアに刺激されました。  
また、テンソルによるテキストの表現については、先行がありましたが、これとは若干コンセプトが異なりました。(R.B. ieee 2004)

# 状態
目次の最適化を記述中です。また、実装のチェックとテストをしています。  
また、衝突検出を準備しています。

# デモ
https://services.limpid-intensity.info/puts.php にサンプルがあります。
擬似的なログインですので、登録は必要ありません。

# ライブラリとしての使い方
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
    std::cout << cstat0.toc(details, detailwords, details, tocwords, tocs, 10);
    
# Tips
これらのプログラムはまだよく検証されたアルゴリズムを使用して*いません*。
abbrev 関数は、理想的には withDetail 関数の逆の最適化となっていることが望ましいですが、現状そうなっていません。  
serialize 関数は、理想的には corpushl(corpus(...)) イニシャライザの逆の最適化となっていることが望ましいですが、現状そうなっていません。  
純粋な興味で作られたプログラムですので、文脈を理解する際の補助として以外には使用しないでください。  
このプログラムを形態要素解析と一緒に使うときちんとした結果が得られますが、現状そうして*いません*。このプログラムは単語同士の距離と順序のみを問題としています、これは何が着目されているかに大まかにつながるものです。  
このプログラムは何が言われたかについて着目するものであり、何が言われていないか、について着目するものではありません。  
もし大きなよく記述されたテキストと辞書があれば、自分の主張と反対の主張を演繹することができますが、これは何が文脈から言えるかであり、何が実際の問題かということではありません。

# その他のダウンロードサイト
* https://ja.osdn.net/projects/puts-cc/
* https://www.sourceforge.net/projects/puts-cc/
* https://konbu.sakura.ne.jp/files/
* https://files.limpid-intensity.info/
