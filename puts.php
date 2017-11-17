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
    exec("mkdir " . $pathb . "output");
    if(!file_exists($pathb . "words.txt"))
      exec("cp ./words.txt " . $pathb);
    file_put_contents($pathb . "email.txt", $_REQUEST["email"]);
  } else
    $pathb = "./datas/" . hash("sha256", $_SESSION["email"] . $_SESSION["salt"]) . "/";
?>
<html>
<head>
<title>Puts sample.</title>
<script language="javascript">
// from: qiita.com/katsunory/items/9bf9ee49ee5c08bf2b3d
function asyncPost(addr, data) {
  var req = new XMLHttpRequest();
  req.onreadystatechange = function () {
    var result = document.getElementById('result');
    if(req.readyState == 4) {
      if(req.status == 200) {
        result.innerHTML = req.responseText + " " + req.status;
      } else {
        result.innerHTML = "Some error occured." + req.status;
      }
    } else {
      result.innerHTML = "Connecting...";
    }
  };
  req.open('POST', addr, true);
  req.setRequestHeader('Content-type',
    'application/x-www-form-urlencoded;charset=UTF-8');
  req.send(data);
  return;
}

function uploadAnalyse() {
  // asyncPost("./analyse.php", "cmd=aa&object=" + encodeURIComponent(document.getElementById("object").value));
  asyncPost("./apply.php", "cmd=aa&containt=" + encodeURIComponent(document.getElementById("object").value));
  return;
}

function updateWords() {
  asyncPost("./apply.php", "cmd=aw&containt=" + encodeURIComponent(document.getElementById("words").value));
  return;
}

function appendTopic() {
  var dom = document.getElementById("new-topic");
  if(document.getElementById("topic-li-" + dom.value) != null) {
    alert("there's a same definition.");
    return;
  }
  asyncPost("./apply.php", "cmd=at&name=" + encodeURIComponent(dom.value) + "&containt=");
  var par   = document.getElementById("topic-ul");
  var last  = document.getElementById("topic-last");
  var newli = document.createElement("li");
  newli.innerHTML = dom.value + " : <br/><textarea maxlength=\"80000\" rows=\"6\" cols=\"40\" id=\"topic-" + dom.value + "\"></textarea><br/><a href=\"javascript: ;\" onClick=\"applyTopic('" + dom.value + "', document.getElementById('topic-" + dom.value + "').value);\">Update</a>|<a href=\"javascript: ;\" onClick=\"deleteTopic('" + dom.value + "');\">-</a>";
  newli.id = "topic-li-" + dom.value;
  par.insertBefore(newli, last);
  dom.value  = "";
  return;
}

function applyTopic(name, containt) {
  asyncPost("./apply.php", "cmd=at&name=" + encodeURIComponent(name) + "&containt=" + encodeURIComponent(containt));
  return;
}

function deleteTopic(name) {
  asyncPost("./apply.php", "cmd=dt&name=" + encodeURIComponent(name));
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
  asyncPost("./apply.php", "cmd=ad&name=" + encodeURIComponent(dom.value) + "&containt=");
  var par   = document.getElementById("dicts-ul");
  var last  = document.getElementById("dicts-last");
  var newli = document.createElement("li");
  newli.innerHTML = dom.value + " : <br/><textarea maxlength=\"80000\" rows=\"6\" cols=\"40\" id=\"dicts-" + dom.value + "\"></textarea><br/><a href=\"javascript: ;\" onClick=\"applyDict('" + dom.value + "', document.getElementById('dicts-" + dom.value + "').value);\">Update</a>|<a href=\"javascript: ;\" onClick=\"deleteDict('" + dom.value + "');\">-</a>";
  newli.id = "dicts-li-" + dom.value;
  par.insertBefore(newli, last);
  dom.value  = "";
  return;
}

function applyDict(name, containt) {
  asyncPost("./apply.php", "cmd=ad&name=" + encodeURIComponent(name) + "&containt=" + encodeURIComponent(containt));
  return;
}

function deleteDict(name) {
  asyncPost("./apply.php", "cmd=dd&name=" + encodeURIComponent(name));
  var dom = document.getElementById("dicts-li-" + name);
  dom.parentElement.removeChild(dom);
  return;
}

function deleteCache() {
  asyncPost("./apply.php", "cmd=da");
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
<?php
  } else {
?>
<p><form action="./puts.php" method="GET"><input type="submit" value="logout" /><input type="hidden" name="logout" /></form>
<a href="javascript:;" onClick="deleteCache();">Delete cache</a>
Your root path: <a href="<?php echo $pathb; ?>">here</a></p>
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
<div style="display:inline-block;vertical-align:top;">
Topics:<br/>
<ul id="topic-ul">
<?php
  foreach (new DirectoryIterator($pathb . 'topics') as $fileInfo) {
    if($fileInfo->isDot()) continue;
    $name = $fileInfo->getFilename();
    echo "<li id=\"topic-li-" . $name . "\">" . $name . " : <br/><textarea maxlength=\"80000\" rows=\"6\" cols=\"40\" id=\"topic-" . $name . "\">";
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
    echo "<li id=\"dicts-li-" . $name . "\">" . $name . " : <br/><textarea maxlength=\"80000\" rows=\"6\" cols=\"40\" id=\"dicts-" . $name . "\">";
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

