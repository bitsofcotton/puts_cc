<?php
function doexecp($so, $se, $cmd, $text, $cwd, $env) {
  $descriptorspec = array(
    0 => array("pipe", "r"),  // stdin.
    1 => array("file", $so, "w"),
    2 => array("file", $se, "w") );
  $process = proc_open($cmd, $descriptorspec, $pipes, $cwd, $env);
  if (is_resource($process)) {
    fwrite($pipes[0], $text . "\n");
    fclose($pipes[0]);
    $return_value = proc_close($process);
  }
  return;
}

$bpath = "./datas/";
foreach(new DirectoryIterator($bpath) as $fileInfo) {
  if(!is_dir($bpath . $fileInfo->getFilename()))
    continue;
  $pathb = $bpath . $fileInfo->getFilename() . "/output/";
  if(!file_exists($pathb))
    continue;
  foreach(new DirectoryIterator($pathb) as $fileInfoSub) {
    $pathb2 = $pathb . $fileInfoSub->getFilename() . "/";
    if($fileInfoSub->isDot())
      continue;
    if(file_exists($pathb2 . "orig.txt") && !file_exists($pathb2 . ".lock")) {
      $cwd = $pathb . '../';
      system('touch ' . $pathb2 . '/.lock');
      $text = "";
      $file = fopen($pathb2 . "orig.txt", "r");
      while(($buf = fgets($file)) !== false)
        $text .= $buf;
      $text .= "\n";
      fclose($file);
      
      doexecp($pathb2 . "lword.txt", $pathb2 . "lword-error.txt",
              "../../puts lword words.txt", $text, $cwd, $env);
      doexecp($pathb2 . "lbalance.txt", $pathb2 . "lbalance-error.txt",
              "../../puts lbalance words.txt", $text, $cwd, $env);
      
      $buf = "";
      foreach (new DirectoryIterator($pathb . '../dicts') as $fileInfo) {
        if($fileInfo->isDot()) continue;
        $name = $fileInfo->getFilename();
        $buf .= "\"dicts/" . $name . "\" ";
      }
      $buf .= " -toc ";
      foreach (new DirectoryIterator($pathb . '../topics') as $fileInfo) {
        if($fileInfo->isDot()) continue;
        $name = $fileInfo->getFilename();
        $buf .= "\"topics/" . $name . "\" ";
      }
      doexecp($pathb2 . "detail.html", $pathb2 . "detail-error.txt",
              "../../puts toc words.txt " . $buf, $text, $cwd, $env);

      $buf = "";
      foreach (new DirectoryIterator($pathb . '../dicts') as $fileInfo) {
        if($fileInfo->isDot()) continue;
        $name = $fileInfo->getFilename();
        $buf .= "\"dicts/" . $name . "\" ";
      }
      doexecp($pathb2 . "stat.html", $pathb2 . "stat-error.txt",
              "../../puts stat words.txt " . $buf, $text, $cwd, $env);
      
      $differs = "";
      $file = fopen($pathb . "../sentry.txt", "r");
      while(($buf = fgets($file)) !== false)
        $differs .= $buf;
      fclose($file);
      if(!preg_match_all("/([0-9a-f]{64,64})/m", $differs, $match)) {
        continue;
      }
      foreach($match[0] as $df) {
        $pathc = $pathb . "../../" . $df . "/";
        if(file_exists($pathc . "dicts")) {
          $buf = " -dict ";
          foreach (new DirectoryIterator($pathb . '../dicts') as $fileInfo) {
            if($fileInfo->isDot()) continue;
            $name = $fileInfo->getFilename();
            $buf .= "\"dicts/" . $name . "\" ";
          }
          $buf .= " -dict2 ";
          foreach (new DirectoryIterator($pathc . '/dicts') as $fileInfo) {
            if($fileInfo->isDot()) continue;
            $name = $fileInfo->getFilename();
            $buf .= "\"" . $pathc . "/dicts/" . $name . "\" ";
          }
          doexecp($pathb2 . $df . "diff.html", $pathb2 . $df . "diff-error.txt",
              "../../puts diff words.txt " . $buf, $text, $cwd, $env);
        }
      }
    }
  }
}
?>
