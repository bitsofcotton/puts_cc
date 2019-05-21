<?php
include("abase.php");

$bpath = "./datas/";
foreach(new DirectoryIterator($bpath) as $fileInfo) {
  if(!is_dir($bpath . $fileInfo->getFilename()) || $fileInfo->isDot())
    continue;
  $pathb = $bpath . $fileInfo->getFilename() . "/output/";
  if(!file_exists($pathb))
    continue;
  foreach(new DirectoryIterator($pathb) as $fileInfoSub) {
    $pathb2 = $pathb . $fileInfoSub->getFilename() . "/";
    if($fileInfoSub->isDot())
      continue;
    if(file_exists($pathb2 . "orig.txt") && !file_exists($pathb2 . ".lock")) {
      system('touch ' . escapeshellarg($pathb2 . '.lock'));
      analyse($pathb . '../', $pathb2 . "orig.txt", true);
    }
  }
}
?>
