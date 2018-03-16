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

function analyse($cwd, $pathb2, $fileb2, $do_stat) {
  echo $pathb2;
  $text = "";
  $file = fopen($fileb2, "r");
  while(($buf = fgets($file)) !== false)
    $text .= $buf;
  fclose($file);
  
  $buf = "";
  foreach (new DirectoryIterator($cwd . '/dicts') as $fileInfo) {
    if($fileInfo->isDot()) continue;
    $name = $fileInfo->getFilename();
    $buf .= "\"dicts/" . $name . "\" ";
  }
  $buf .= " -toc ";
  foreach (new DirectoryIterator($cwd . '/topics') as $fileInfo) {
    if($fileInfo->isDot()) continue;
    $name = $fileInfo->getFilename();
    $buf .= "\"topics/" . $name . "\" ";
  }
  doexecp($pathb2 . "-detail.html", $pathb2 . "-detail-error.txt",
          "../../puts toc words.txt " . $buf, $text, $cwd, $env);
  
  if($do_stat) {
    $buf = "";
    foreach (new DirectoryIterator($cwd . '/dicts') as $fileInfo) {
      if($fileInfo->isDot()) continue;
      $name = $fileInfo->getFilename();
      $buf .= "\"dicts/" . $name . "\" ";
    }
    doexecp($pathb2 . "-stat.html", $pathb2 . "-stat-error.txt",
            "../../puts stat words.txt " . $buf, $text, $cwd, $env);
  }
  doexecp($pathb2 . "-lbalance.txt", $pathb2 . "-lbalance-error.txt",
          "../../puts lbalance words.txt", $text, $cwd, $env);
  doexecp($pathb2 . "-lword.txt", $pathb2 . "-lword-error.txt",
          "../../puts lword words.txt", $text, $cwd, $env);
  
  $differs = "";
  $file = fopen($cwd . "/sentry.txt", "r");
  while(($buf = fgets($file)) !== false)
    $differs .= $buf;
  fclose($file);
  if(!preg_match_all("/([0-9a-f]{64,64})/m", $differs, $match)) {
    return;
  }
  foreach($match[0] as $df) {
    $pathc = $cwd . "../" . $df . "/";
    if(file_exists($pathc . "dicts")) {
      $buf = " -dict ";
      foreach (new DirectoryIterator($cwd . '/dicts') as $fileInfo) {
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
  return;
}
?>
