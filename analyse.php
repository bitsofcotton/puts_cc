<?php
foreach(new DirectoryIterator('./datas/') as $fileInfo) {
  if(!is_dir('./datas/' . $fileInfo->getFilename()))
    continue;
  $pathb = "./datas/" . $fileInfo->getFilename() . "/output/";
  if(!file_exists($pathb))
    continue;
  foreach(new DirectoryIterator($pathb) as $fileInfoSub) {
    if($fileInfoSub->isDot())
      continue;
    $pathb2 = $pathb . $fileInfoSub->getFilename() . '/';
    $cwd    = $pathb . '../';
    if(file_exists($pathb2 . "orig-ref.txt") && !file_exists($pathb2 . ".lock")) {
      echo $pathb2;
      system('touch ' . $pathb2 . '/.lock');
      $name = $pathb2 . "orig-ref.txt";
      $f = fopen($pathb2 . "orig-ref-dicts.txt", "r");
      $pathc = fgets($f);
      fclose($f);
      $text = "";
      $file = fopen($name, "r");
      while(($buf = fgets($file)) !== false)
        $text .= $buf;
      fclose($file);
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
        $buf .= "\"../../" . $pathc . "dicts/" . $name . "\" ";
      }
      $descriptorspec = array(
        0 => array("pipe", "r"),  // stdin.
        1 => array("file", $pathb2 . "diff.html", "w"),  // stdout.
        2 => array("file", $pathb2 . "diff-error.txt", "w") // stderr.
      );
      $f2 = fopen($pathb2 . "comm.txt", "w");
      fwrite($f2, "cmd:" . $buf);
      fclose($f2);
      $process = proc_open('../../puts diff words.txt ' . $buf, $descriptorspec, $pipes, $cwd, $env);
      if (is_resource($process)) {
        fwrite($pipes[0], $text . "\n");
        fclose($pipes[0]);
        $return_value = proc_close($process);
      }
    } else if(file_exists($pathb2 . "orig.txt") && !file_exists($pathb2 . ".lock")) {
      system('touch ' . $pathb2 . '/.lock');
      $name = $pathb2 . "orig.txt";
      $text = "";
      $file = fopen($name, "r");
      while(($buf = fgets($file)) !== false)
        $text .= $buf;
      fclose($file);
      $descriptorspec = array(
        0 => array("pipe", "r"),  // stdin.
        1 => array("file", $pathb2 . "lword.txt", "w"),  // stdout.
        2 => array("file", $pathb2 . "lword-error.txt", "w") // stderr.
      );
      $process = proc_open('../../puts lword', $descriptorspec, $pipes, $cwd, $env);
      if (is_resource($process)) {
        fwrite($pipes[0], $text . "\n");
        fclose($pipes[0]);
        $return_value = proc_close($process);
      }
      
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
      
      $descriptorspec = array(
        0 => array("pipe", "r"),  // stdin.
        1 => array("file", $pathb2 . "detail.txt", "w"),  // stdout.
        2 => array("file", $pathb2 . "detail-error.txt", "w") // stderr.
      );
      $process = proc_open('../../puts toc words.txt ' . $buf, $descriptorspec, $pipes, $cwd, $env);
      
      if (is_resource($process)) {
        fwrite($pipes[0], $text . "\n");
        fclose($pipes[0]);
        $return_value = proc_close($process);
      }

      $buf = "";
      foreach (new DirectoryIterator($pathb . '../dicts') as $fileInfo) {
        if($fileInfo->isDot()) continue;
        $name = $fileInfo->getFilename();
        $buf .= "\"dicts/" . $name . "\" ";
      }
      
      $descriptorspec = array(
        0 => array("pipe", "r"),  // stdin.
        1 => array("file", $pathb2 . "stat.html", "w"),  // stdout.
        2 => array("file", $pathb2 . "stat-error.txt", "w") // stderr.
      );
      $process = proc_open('../../puts stat words.txt ' . $buf, $descriptorspec, $pipes, $cwd, $env);
      if (is_resource($process)) {
        fwrite($pipes[0], $text . "\n");
        fclose($pipes[0]);
        $return_value = proc_close($process);
      }
      
      $descriptorspec = array(
        0 => array("pipe", "r"),  // stdin.
        1 => array("file", $pathb2 . "redig.txt", "w"),  // stdout.
        2 => array("file", $pathb2 . "redig-error.txt", "w") // stderr.
      );
      $process = proc_open('../../puts redig words.txt ', $descriptorspec, $pipes, $cwd, $env);
      if (is_resource($process)) {
        fwrite($pipes[0], $text . "\n");
        fclose($pipes[0]);
        $return_value = proc_close($process);
      }
    }
  }
}

foreach(new DirectoryIterator('./datas/') as $fileInfo) {
  if(!is_dir('./datas/' . $fileInfo->getFilename()))
    continue;
  $pathb = "./datas/" . $fileInfo->getFilename() . "/web/";
  if(!file_exists($pathb))
    continue;
  foreach(new DirectoryIterator($pathb) as $fileInfoSub) {
    if($fileInfoSub->isDot())
      continue;
    $pathb2 = $pathb . $fileInfoSub->getFilename() . '/';
    $cwd    = $pathb . '../';
    preg_match("/\d+\.txt$/", $pathb2, $match);
    if(count($match) > 0 && !file_exists($pathb2 . "-redig.txt")) {
      echo $pathb2;
      $text = "";
      $file = fopen($pathb2, "r");
      while(($buf = fgets($file)) !== false)
        $text .= $buf;
      fclose($file);
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
      
      $descriptorspec = array(
        0 => array("pipe", "r"),  // stdin.
        1 => array("file", $pathb2 . "-summary.txt", "w"),  // stdout.
        2 => array("file", $pathb2 . "-summary-error.txt", "w") // stderr.
      );
      $process = proc_open('../../puts toc words.txt ' . $buf, $descriptorspec, $pipes, $cwd, $env);
      if (is_resource($process)) {
        fwrite($pipes[0], $text . "\n");
        fclose($pipes[0]);
        $return_value = proc_close($process);
      }
      
      $buf = "";
      foreach (new DirectoryIterator($pathb . '../dicts') as $fileInfo) {
        if($fileInfo->isDot()) continue;
        $name = $fileInfo->getFilename();
        $buf .= "\"dicts/" . $name . "\" ";
      }
      
      $descriptorspec = array(
        0 => array("pipe", "r"),  // stdin.
        1 => array("file", $pathb2 . "-stat.txt", "w"),  // stdout.
        2 => array("file", $pathb2 . "-stat-error.txt", "w") // stderr.
      );
      $process = proc_open('../../puts stat words.txt ' . $buf, $descriptorspec, $pipes, $cwd, $env);
      if (is_resource($process)) {
        fwrite($pipes[0], $text . "\n");
        fclose($pipes[0]);
        $return_value = proc_close($process);
      }
      
      $descriptorspec = array(
        0 => array("pipe", "r"),  // stdin.
        1 => array("file", $pathb2 . "-redig.txt", "w"),  // stdout.
        2 => array("file", $pathb2 . "-redig-error.txt", "w") // stderr.
      );
      $process = proc_open('../../puts redig words.txt', $descriptorspec, $pipes, $cwd, $env);
      if (is_resource($process)) {
        fwrite($pipes[0], $text . "\n");
        fclose($pipes[0]);
        $return_value = proc_close($process);
      }
    }
  }
}
?>
