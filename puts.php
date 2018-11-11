<?php
  session_start();
  if(isset($_REQUEST["logout"])) {
    session_destroy();
    header("location: puts.php");
    exit();
  }
  if(isset($_REQUEST["email"]) || isset($_REQUEST["salt"])) {
    if(isset($_REQUEST["email"]))
      $_SESSION["email"] = $_REQUEST["email"];
    if(isset($_REQUEST["salt"]))
      $_SESSION["salt"]  = $_REQUEST["salt"];
    $pathb = "./datas/" . hash("sha256", $_SESSION["email"] . $_SESSION["salt"]) . "/";
    exec("mkdir " . $pathb);
    exec("mkdir " . $pathb . "topics");
    exec("mkdir " . $pathb . "dicts");
    exec("mkdir " . $pathb . "pdict");
    exec("mkdir " . $pathb . "output");
    exec("mkdir " . $pathb . "web");
    exec("mkdir " . $pathb . "web/sub");
    exec("mkdir " . $pathb . "crawl");
    exec("rm -f " . $pathb . "puts.core");
    if(!file_exists($pathb . "words.txt"))
      exec("cp ./words.txt " . $pathb);
    exec("cp ./style.css " . $pathb);
    if(!file_exists($pathb . "sentry.txt"))
      exec("touch " . $pathb . "sentry.txt");
    exec("touch " . $pathb . "urls.txt");
    file_put_contents($pathb . "email.txt", $_REQUEST["email"]);
  } else
    $pathb = "./datas/" . hash("sha256", $_SESSION["email"] . $_SESSION["salt"]) . "/";
?>
<html>
<head>
<title>Puts sample.</title>
<script language="javascript">
// thanks: qiita.com/katsunory/items/9bf9ee49ee5c08bf2b3d
function asyncPost(addr, data) {
  var req = new XMLHttpRequest();
  req.onreadystatechange = function () {
    var result = document.getElementById('result');
    if(req.readyState == 4) {
      if(req.status == 200) {
        result.innerHTML = req.responseText + " " + req.status;
      } else {
        result.innerHTML = "Some error had be occur." + req.status;
      }
    } else {
      result.innerHTML = "Connecting...";
    }
  };
  req.open('POST', addr, true);
  req.send(data);
  return;
}

function uploadAnalyse() {
  fd = new FormData();
  fd.append("cmd", "aa");
  fd.append("containt", document.getElementById("object").value);
  asyncPost("./apply.php", fd);
  return;
}

function updateWords() {
  fd = new FormData();
  fd.append("cmd", "aw");
  fd.append("containt", document.getElementById("words").value);
  asyncPost("./apply.php", fd);
  return;
}

function updateSentryList() {
  fd = new FormData();
  fd.append("cmd", "st");
  fd.append("containt", document.getElementById("sentry_list").value);
  asyncPost("./apply.php", fd);
  return;
}

function updateURLs() {
  fd = new FormData();
  fd.append("cmd", "url");
  fd.append("containt", document.getElementById("urls").value);
  asyncPost("./apply.php", fd);
  return;
}

function appendTopic() {
  var dom = document.getElementById("new-topic");
  if(document.getElementById("topic-li-" + dom.value) != null) {
    alert("there's a same definition.");
    return;
  }
  fd = new FormData();
  fd.append("cmd", "at");
  fd.append("name", dom.value);
  fd.append("containt", "");
  asyncPost("./apply.php", fd);
  var par   = document.getElementById("topic-ul");
  var last  = document.getElementById("topic-last");
  var newli = document.createElement("li");
  newli.innerHTML = dom.value + " : <br/><textarea maxlength=\"80000\" rows=\"12\" cols=\"80\" id=\"topic-" + dom.value + "\"></textarea><br/><a href=\"javascript: ;\" onClick=\"applyTopic('" + dom.value + "', document.getElementById('topic-" + dom.value + "').value);\">Update</a>|<a href=\"javascript: ;\" onClick=\"deleteTopic('" + dom.value + "');\">-</a>";
  newli.id = "topic-li-" + dom.value;
  par.insertBefore(newli, last);
  dom.value  = "";
  return;
}

function applyTopic(name, containt) {
  fd = new FormData();
  fd.append("cmd", "at");
  fd.append("name", name);
  fd.append("containt", containt);
  asyncPost("./apply.php", fd);
  return;
}

function deleteTopic(name) {
  fd = new FormData();
  fd.append("cmd", "dt");
  fd.append("name", name);
  asyncPost("./apply.php", fd);
  var dom = document.getElementById("topic-li-" + name);
  dom.parentElement.removeChild(dom);
  return;
}

function appendDict() {
  var dom = document.getElementById("new-dicts");
  if(document.getElementById("dicts-li-" + dom.value) != null) {
    alert("there's a same definition.");
    return;
  }
  fd = new FormData();
  fd.append("cmd", "ad");
  fd.append("name", dom.value);
  fd.append("containt", "");
  asyncPost("./apply.php", fd);
  var par   = document.getElementById("dicts-ul");
  var last  = document.getElementById("dicts-last");
  var newli = document.createElement("li");
  newli.innerHTML = dom.value + " : <br/><textarea maxlength=\"80000\" rows=\"12\" cols=\"80\" id=\"dicts-" + dom.value + "\"></textarea><br/><a href=\"javascript: ;\" onClick=\"applyDict('" + dom.value + "', document.getElementById('dicts-" + dom.value + "').value);\">Update</a>|<a href=\"javascript: ;\" onClick=\"deleteDict('" + dom.value + "');\">-</a>";
  newli.id = "dicts-li-" + dom.value;
  par.insertBefore(newli, last);
  dom.value  = "";
  return;
}

function applyDict(name, containt) {
  fd = new FormData();
  fd.append("cmd", "ad");
  fd.append("name", name);
  fd.append("containt", containt);
  asyncPost("./apply.php", fd);
  return;
}

function deleteDict(name) {
  fd = new FormData();
  fd.append("cmd", "dd");
  fd.append("name", name);
  asyncPost("./apply.php", fd);
  var dom = document.getElementById("dicts-li-" + name);
  dom.parentElement.removeChild(dom);
  return;
}

function deleteCache() {
  fd = new FormData();
  fd.append("cmd", "da");
  asyncPost("./apply.php", fd);
  return;
}

function archive() {
  fd = new FormData();
  fd.append("cmd", "ta");
  asyncPost("./apply.php", fd);
  return;
}
</script>
</head>
<body>
<div id="result"></div>
<div>
<?php
  if(!isset($_SESSION["email"]) || !isset($_SESSION["salt"])) {
?>
<p>
  <form action="./puts.php" method="GET">
  <ul>
  <li>Email: <input type="email" name="email" /></li>
  <li>Salt:  <input type="text" name="salt" /></li>
  <li><input type="submit" value="Login" /></li>
  </ul>
  </form>
</p>
<?php
  } else if(isset($_REQUEST["adddict"])) {
?>
<p>
  Add dictionary Edit:
<?php
  {
    $name = $_REQUEST["name"];
    echo "<ul id=\"dicts-ul\">";
    echo "<li>";
    echo "Please edit and paste:";
    echo "<textarea maxlength=\"80000\" rows=\"30\" cols=\"90\" id=\"dicts-" . $name . "\">";
    echo $_REQUEST["entry"] . "\n";
    echo "</textarea></li>";
    echo "<li id=\"dicts-last\"><input type=\"text\" id=\"new-dicts\" /><a href=\"javascript: ;\" onClick=\"appendDict();\">+</a></li>";
    echo "</ul>";
  }
?>
</p>
<?php
  } else {
?>
<p>
<form action="./puts.php" method="GET"><input type="submit" value="logout" /><input type="hidden" name="logout" /></form>
<a href="javascript:;" onClick="deleteCache();">Delete cache</a> |
<a href="<?php echo $pathb; ?>">Your root path</a> |
<a href="<?php echo $pathb; ?>/web">Per 6 hours news</a> |
<a href="javascript:;" onClick="archive();">Archive all</a></p>
</p>
<p>
Analyse text:<br/>
<textarea maxlength="80000" rows="12" cols="80" name="object" id="object"></textarea><br/>
<a href="javascript:;" onClick="uploadAnalyse();">upload analyse base.</a>
</p>
<p>
Word list:<br/>
<textarea rows="12" cols="80" name="words" id="words">
<?php
  $file = fopen($pathb . "words.txt", "r");
  while(($buf = fgets($file)) !== false) {
    echo $buf;
  }
  fclose($file);
?></textarea><br/>
<a href="javascript: ;" onClick="updateWords();">update</a>
</p>
<p>
Scrap URL list:<br/>
<textarea rows="12" cols="80" name="urls" id="urls">
<?php
  $file = fopen($pathb . "urls.txt", "r");
  while(($buf = fgets($file)) !== false) {
    echo $buf;
  }
  fclose($file);
?></textarea><br/>
<a href="javascript: ;" onClick="updateURLs();">update</a>
</p>
<p>
Sentry Dictionary Directory List (Hash only):<br/>
<textarea rows="8" cols="80" name="sentry_list" id="sentry_list">
<?php
  $file = fopen($pathb . "sentry.txt", "r");
  while(($buf = fgets($file)) !== false) {
    echo $buf;
  }
  fclose($file);
?></textarea><br/>
<a href="javascript: ;" onClick="updateSentryList();">update</a>
</p>
<div style="display:inline-block;vertical-align:top;">
Topics:<br/>
<ul id="topic-ul">
<?php
  foreach (new DirectoryIterator($pathb . 'topics') as $fileInfo) {
    if($fileInfo->isDot()) continue;
    $name = $fileInfo->getFilename();
    echo "<li id=\"topic-li-" . $name . "\">" . $name . " : <br/><textarea maxlength=\"80000\" rows=\"12\" cols=\"80\" id=\"topic-" . $name . "\">";
    $file = fopen($pathb . "topics/" . $name, "r");
    while(($buf = fgets($file)) !== false) {
      echo $buf;
    }
    fclose($file);
    echo "</textarea><br/>";
    echo "<a href=\"javascript: ;\" onClick=\"applyTopic('" . $name . "', document.getElementById('topic-" . $name . "').value);\">Update</a>|";
    echo "<a href=\"javascript: ;\" onClick=\"deleteTopic('" . $name . "');\">-</a></li>";
  }
  echo "<li id=\"topic-last\"><input type=\"text\" id=\"new-topic\" /><a href=\"javascript: ;\" onClick=\"appendTopic();\">+</a></li>";
?>
</ul>
</div>
<div style="display:inline-block;vertical-align:top;">
Dictionaries:<br/>
<ul id="dicts-ul">
<?php
  foreach (new DirectoryIterator($pathb . 'dicts') as $fileInfo) {
    if($fileInfo->isDot()) continue;
    $name = $fileInfo->getFilename();
    echo "<li id=\"dicts-li-" . $name . "\">" . $name . " : <br/><textarea maxlength=\"80000\" rows=\"12\" cols=\"80\" id=\"dicts-" . $name . "\">";
    $file = fopen($pathb . "dicts/" . $name, "r");
    while(($buf = fgets($file)) !== false) {
      echo $buf;
    }
    fclose($file);
    echo "</textarea><br/>";
    echo "<a href=\"javascript: ;\" onClick=\"applyDict('" . $name . "', document.getElementById('dicts-" . $name . "').value);\">Update</a>|";
    echo "<a href=\"javascript: ;\" onClick=\"deleteDict('" . $name . "');\">-</a></li>";
  }
  echo "<li id=\"dicts-last\"><input type=\"text\" id=\"new-dicts\" /><a href=\"javascript: ;\" onClick=\"appendDict();\">+</a></li>";
?>
</ul>
</div>
<?php
  }
?>

</div>
</body>
</html>

