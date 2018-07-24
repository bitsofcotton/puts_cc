<?php
session_start();
$M_COUNT = 300;
$M_STR   = 80000;
if(!isset($_SESSION["email"]) || !isset($_SESSION["salt"])) {
  echo "Please login.";
  return;
}
$pathb = "./datas/" . hash("sha256", $_SESSION["email"] . $_SESSION["salt"]) . "/";
$name  = "";
if(isset($_REQUEST["name"]))
  $name = basename($_REQUEST["name"]);
$containt = "";
if(isset($_REQUEST["containt"]))
  $containt = $_REQUEST["containt"];
if(!isset($_REQUEST["cmd"])) {
  echo "illegal command.";
  return;
}
$cmd = $_REQUEST["cmd"];
switch($cmd) {
case "aa":
  if(mb_strlen($containt) < $M_STR) {
    $pathb2 = 'output/' . hash("sha256", $containt) . "/";
    exec("mkdir -p " . $pathb . $pathb2);
    file_put_contents($pathb . $pathb2 . "/orig.txt", $containt);
    exec("rm -f " . $pathb . $pathb2 . "/.lock");
    echo "Your path is <a href=\"" . $pathb . $pathb2 . "\">Here</a>";
    exec("php ./analyse.php");
  }
  break;
case "aw":
  if($containt == "") {
    exec("rm -f " . $pathb . "words.txt");
  } else {
    file_put_contents($pathb . "words.txt", $containt);
  }
  break;
case "url":
  file_put_contents($pathb . "urls.txt", $containt);
  break;
case "st":
  $pathc = "";
  $ctr   = 0;
  if(preg_match_all("/([0-9a-f]{64,64})/m", $containt, $match)) {
    foreach($match[0] as $m) {
      if(file_exists($pathb . "../" . $m)) {
        $pathc .= $m . "\n";
        $ctr   += 1;
        if($ctr >= 8) {
           break;
        }
      }
    }
  }
  file_put_contents($pathb . "sentry.txt", $pathc);
  break;
case "at":
  $fi = new FilesystemIterator($pathb . "topics/", FilesystemIterator::SKIP_DOTS);
  if(iterator_count($fi) < $M_COUNT && mb_strlen($containt) < $M_STR)
    file_put_contents($pathb . "topics/" . $name, $containt);
  break;
case "dt":
  unlink($pathb . "topics/" . $name);
  break;
case "ad":
  $fi = new FilesystemIterator($pathb . "dicts/", FilesystemIterator::SKIP_DOTS);
  if(iterator_count($fi) < $M_COUNT && mb_strlen($containt) < $M_STR)
    file_put_contents($pathb . "dicts/" . $name, $containt);
  break;
case "dd":
  unlink($pathb . "dicts/" . $name);
  break;
case "da":
  exec("rm -fr " . $pathb . "/output && mkdir -p " . $pathb . "/output");
  break;
case "ta":
  exec("cd " . $pathb . " && rm -f *.core");
  exec("cd " . $pathb . " && rm -f archive.tar.gz && tar cfz archive.tar.gz .");
  break;
default:
  echo "invalid cmd.";
  return;
}
?>
