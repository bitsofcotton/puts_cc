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

function prepdicts($cwd, $txt, $text, $env) {
  doexecp($txt . "-prep.txt", $txt . "-prep-err.txt",
          "../../puts prep words.txt", $text, $cwd, $env);
  $f  = fopen($txt . "-prep.txt", "r");
  while(($buf = fgets($f)) !== false) {
    doexecp($cwd . "/pdict/" . basename(chop($buf)), "/dev/null",
            "python ../../prep.py " . escapeshellarg(basename(chop($buf))),
            "\n", $cwd, $env);
  }
  doexecp("/dev/null", "/dev/null", "sh -c 'find pdict/ -empty | xargs rm'", "\n", $cwd, $env);
  fclose($f);
}

function analyse($cwd, $pathb2, $do_stat) {
  echo $pathb2;
  $text = "";
  $file = fopen($pathb2, "r");
  while(($buf = fgets($file)) !== false)
    $text .= $buf;
  fclose($file);
  
  doexecp("/dev/null", "/dev/null", "rm -f pdict/*",
          "\n", $cwd, $env);
  prepdicts($cwd, $pathb2, $text, $env);
  foreach (new DirectoryIterator($cwd . '/topics') as $fileInfo) {
    if($fileInfo->isDot()) continue;
    $ltext = "";
    $lf = fopen($cwd . '/topics/' . $fileInfo->getFilename(), "r");
    while(($buf = fgets($lf)) !== false)
      $ltext .= $buf;
    fclose($lf);
    prepdicts($cwd, $pathb2, $ltext, $env);
  }
  
  $buf = "";
  foreach (new DirectoryIterator($cwd . '/dicts') as $fileInfo) {
    if($fileInfo->isDot()) continue;
    $name = $fileInfo->getFilename();
    $buf .= escapeshellarg("dicts/" . $name) . " ";
  }
  foreach (new DirectoryIterator($cwd . '/pdict') as $fileInfo) {
    if($fileInfo->isDot()) continue;
    $name = $fileInfo->getFilename();
    $buf .= escapeshellarg("pdict/" . $name) . " ";
  }
  $buf .= " -toc ";
  foreach (new DirectoryIterator($cwd . '/topics') as $fileInfo) {
    if($fileInfo->isDot()) continue;
    $name = $fileInfo->getFilename();
    $buf .= escapeshellarg("topics/" . $name) . " ";
  }
  doexecp($pathb2 . "-detail.html", $pathb2 . "-detail-error.txt",
          "../../puts toc words.txt " . $buf, $text, $cwd, $env);
  doexecp($pathb2 . "-lack.html", $pathb2 . "-lack-error.txt",
          "../../puts lack words.txt " . $buf, $text, $cwd, $env);
  
  if($do_stat) {
    $buf = "";
    foreach (new DirectoryIterator($cwd . '/dicts') as $fileInfo) {
      if($fileInfo->isDot()) continue;
      $name = $fileInfo->getFilename();
      $buf .= escapeshellarg("dicts/" . $name) . " ";
    }
    foreach (new DirectoryIterator($cwd . '/pdict') as $fileInfo) {
      if($fileInfo->isDot()) continue;
      $name = $fileInfo->getFilename();
      $buf .= escapeshellarg("pdict/" . $name) . " ";
    }
    doexecp($pathb2 . "-stat.html", $pathb2 . "-stat-error.txt",
            "../../puts stat words.txt " . $buf, $text, $cwd, $env);
    doexecp($pathb2 . "-root.html", $pathb2 . "-root-error.txt",
            "../../puts findroot words.txt " . $buf, $text, $cwd, $env);
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
    $pathc = "../" . $df . "/";
    if(file_exists($cwd . '/' . $pathc . "dicts")) {
      $buf  = " -dict ";
      foreach (new DirectoryIterator($cwd . '/dicts') as $fileInfo) {
        if($fileInfo->isDot()) continue;
        $name = $fileInfo->getFilename();
        $buf .= escapeshellarg("dicts/" . $name) . " ";
      }
      foreach (new DirectoryIterator($cwd . '/pdict') as $fileInfo) {
        if($fileInfo->isDot()) continue;
        $name = $fileInfo->getFilename();
        $buf .= escapeshellarg("pdict/" . $name) . " ";
      }
      $buf .= " -dict2 ";
      foreach (new DirectoryIterator($cwd . '/' . $pathc . '/dicts') as $fileInfo) {
        if($fileInfo->isDot()) continue;
        $name = $fileInfo->getFilename();
        $buf .= escapeshellarg($pathc . "/dicts/" . $name) . " ";
      }
      foreach (new DirectoryIterator($cwd . '/pdict') as $fileInfo) {
        if($fileInfo->isDot()) continue;
        $name = $fileInfo->getFilename();
        $buf .= escapeshellarg("pdict/" . $name) . " ";
      }
      doexecp($pathb2 . $df . "diff.html", $pathb2 . $df . "diff-error.txt",
          "../../puts diff words.txt " . $buf, $text, $cwd, $env);
    }
  }
  return;
}
?>
