<?php

// Usage:
//  move *.php prep.py puts into the directory.
//  mkdir data directory and configure permissions.
//  we might need python-mysql environments.
//
//  prep.py needs legacy _mysql with python2.7, if we need, please rewrite them.
//
//  we need mediawiki and xml2sql softwared to import wiktionary database (GFDL).
//  however dictating wiktionary needs huge amounts of memory and very slow.

$dont_pdict = 1;

$dameji = array("!", '"', "#", "$", "%", "&", "'", "(", ")", "*", "+", ",", "-", "/", ":", ";", "<", "=", ">", "?", "@", "[", "\\", "]", "^", "_", "`", "{", "|", "}", "~", " ");
$damejinodir = array("!", '"', "#", "$", "%", "&", "'", "(", ")", "*", "+", ",", "-", ":", ";", "<", "=", ">", "?", "@", "[", "\\", "]", "^", "_", "`", "{", "|", "}", "~", " ");

function prepdicts($cwd, $text) {
  global $dont_pdict;
  global $dameji;
  if($dont_pdict)
    return;
  exec("./puts prep words.txt < " . $text . " > " . $text . "-prep.txt");
  $f   = fopen($text . "-prep.txt", "r");
  $buf = "";
  while(($buf = fgets($f)) !== false) {
    $buf .= str_replace($dameji, "", $buf) . " ";
  }
  fclose($f);
  exec("python3 prep.py " . $buf);
  exec("find pdict/ -empty | xargs rm -f");
  return;
}

function analyse($cwd, $text) {
  global $dameji;
  echo $text;
  // prepare dicts 
  exec("rm -f pdict/*");
  prepdicts($cwd, $text);
  foreach (new DirectoryIterator($cwd . '/topics') as $fileInfo) {
    if($fileInfo->isDot()) continue;
    $name = $fileInfo->getFilename();
    $lf = $cwd . '/topics/' . str_replace($dameji, "", $name);
    prepdicts($cwd, $lf);
  }
  foreach (new DirectoryIterator($cwd . '/dicts') as $fileInfo) {
    if($fileInfo->isDot()) continue;
    $name = $fileInfo->getFilename();
    $lf = $cwd . '/dicts/' . str_replace($dameji, "", $name);
    prepdicts($cwd, $lf);
  }
  
  $buf = "";
  foreach (new DirectoryIterator($cwd . '/dicts') as $fileInfo) {
    if($fileInfo->isDot()) continue;
    $name = $fileInfo->getFilename();
    $buf .= $cwd . "/dicts/" . str_replace($dameji, "", $name) . " ";
  }
  foreach (new DirectoryIterator($cwd . '/pdict') as $fileInfo) {
    if($fileInfo->isDot()) continue;
    $name = $fileInfo->getFilename();
    $buf .= $cwd . "/pdict/" . str_replace($dameji, "", $name) . " ";
  }
  $buf2 = $buf . " -toc ";
  foreach (new DirectoryIterator($cwd . '/topics') as $fileInfo) {
    if($fileInfo->isDot()) continue;
    $name = $fileInfo->getFilename();
    $buf2 .= $cwd . "/topics/" . str_replace($dameji, "", $name) . " ";
  }
  exec("./puts toc  words.txt " . $buf2 . " < " . $text . " > " . $text . "-detail.html");
  exec("./puts lack words.txt " . $buf2 . " < " . $text . " > " . $text . "-lack.html");
  exec("./puts stat words.txt " . $buf . " < " . $text . " > " . $text . "-stat.html");
  exec("./puts findroot words.txt " . $buf . " < " . $text . " > " . $text . "-root.html");
  exec("./puts lword words.txt < " . $text . " > " . $text . "-lword.txt");
  exec("./puts lbalance words.txt < " . $text . " > " . $text . "-lbalance.txt");
  $differs = "";
  $file = fopen($cwd . "/sentry.txt", "r");
  while(($buf = fgets($file)) !== false)
    $differs .= $buf;
  fclose($file);
  if(!preg_match_all("/([0-9a-f]{64,64})/m", $differs, $match)) {
    return;
  }
  foreach($match[0] as $df) {
    echo $df;
    $pathc = "../" . $df . "/";
    if(file_exists($cwd . '/' . $pathc . "dicts")) {
      $buf3 = " -dict2 ";
      foreach (new DirectoryIterator($cwd . "/" . $pathc . "dicts") as $fileInfo) {
        if($fileInfo->isDot()) continue;
        $name = $fileInfo->getFilename();
        $buf3 .= $cwd . "/" . $pathc . "dicts/" . str_replace($dameji, "", $name) . " ";
      }
      foreach (new DirectoryIterator($cwd . "/" . $pathc . "pdict") as $fileInfo) {
        if($fileInfo->isDot()) continue;
        $name = $fileInfo->getFilename();
        $buf3 .= $cwd . "/" . $pathc . "pdict/" . str_replace($dameji, "", $name) . " ";
      }
      $buf3 = $buf . $buf3;
      exec("./puts diff words.txt " . $buf3 . " < " . $text . " > " . $text . $df . "diff.html");
      exec("./puts same words.txt " . $buf3 . " < " . $text . " > " . $text . $df . "same.html");
    }
  }
  return;
}

$bpath = "/var/www/htdocs/datas/";
foreach(new DirectoryIterator($bpath) as $fileInfo) {
  if(!is_dir($bpath . $fileInfo->getFilename()) || $fileInfo->isDot())
    continue;
  $pathb = $bpath . $fileInfo->getFilename() . "/output/";
  if(file_exists($pathb)) {
    foreach(new DirectoryIterator($pathb) as $fileInfoSub) {
      $pathb2 = $pathb . $fileInfoSub->getFilename() . "/";
      if($fileInfoSub->isDot())
        continue;
      $pathb2 = str_replace($damejinodir, "", $pathb2);
      if(strpos($pathb2, ".lock"))
        continue;
      if(file_exists($pathb2 . "orig.txt") && !file_exists($pathb2 . ".lock")) {
        system('touch ' . $pathb2 . ".lock");
        analyse($pathb . '../', $pathb2 . "orig.txt");
      }
    }
  }
  $pathb = $bpath . $fileInfo->getFilename() . "/crawl/";
  if(file_exists($pathb)) {
    foreach(new DirectoryIterator($pathb) as $fileInfoSub) {
      $pathb2 = $pathb . $fileInfoSub->getFilename();
      if($fileInfoSub->isDot())
        continue;
      if(strpos($pathb2, ".lock"))
        continue;
      $pathb2 = str_replace($damejinodir, "", $pathb2);
      if(strlen(basename($pathb2)) == 14 &&
         file_exists($pathb2) && !file_exists($pathb2 . ".lock")) {
        system('touch ' . $pathb2 . '.lock');
        analyse($pathb . '../',  $pathb2);
      }
    }
  }
}
?>
