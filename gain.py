import mechanize
import re
import ipaddress
import socket

mechanize._sockettimeout._GLOBAL_DEFAULT_TIMEOUT = 100

ipranges = [[u"150.147.0.0", u"150.147.255.255"],
            [u"150.246.0.0", u"150.246.255.255"],
            [u"150.249.0.0", u"150.249.255.255"],
            [u"150.120.0.0", u"150.120.127.255"],
            [u"153.120.128.0", u"153.120.191.255"],
            [u"153.120.192.0", u"153.120.221.255"],
            [u"153.120.222.0", u"153.120.223.255"],
            [u"153.120.224.0", u"153.120.255.255"],
            [u"153.121.0.0",   u"153.121.31.255"],
            [u"153.121.32.0",  u"153.121.63.255"],
            [u"153.121.64.0",  u"152.121.95.255"],
            [u"153.121.96.0",  u"153.121.127.255"]]
outbase  = '~/Downloads'

urls  = re.compile(r'href=[\"\']([^\"\']+)[\"\']')
urlsr = re.compile(r"([a-zA-Z0-9\.\-\_]+(:[0-9]+)?/[a-zA-Z0-9\=\?\~/\-\_\.\&]+)")
replscriptstart = re.compile(r"<script")
replscriptend   = re.compile(r"/script>")
repltags        = re.compile(r"<[^<>]+>")

def initBrowser():
  br = mechanize.Browser()
  br.set_handle_equiv(True)
  br.set_handle_redirect(True)
  br.set_handle_referer(True)
  br.user_agent_alias = 'Windows IE 9'
  br.set_handle_robots(False)
  return br

def cutLinks(html):
  urlf = urls.findall(html)
  urlfr = urlsr.findall(html)
  items = []
  for item in urlf:
    items.append(item)
  for item in urlfr:
    items.append(item[0])
  return items

def innerHTMLs(html):
  scstarts = replscriptstart.finditer(html)
  scends   = replscriptend.finditer(html)
  sanitized = ""
  starts    = 0
  for scs in scstarts:
    scec = 0
    for sce in scends:
      if(scs.start() < sce.end()):
        scec = sce.end()
        break
    if(scec > 0):
      sanitized += html[starts:scs.start()]
      starts = scec
  sanitized += html[starts:- 1]
  repltags = re.compile(r"<[^<>]+>")
  return re.sub(repltags, "", sanitized)

def ping80(host):
  sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
  sock.settimeout(1)
  try:
    result = sock.connect_ex((host, 80))
  except:
    False
  if result != 0:
    return False
  return True

for iprange in ipranges:
  start  = int(ipaddress.IPv4Address(iprange[0]))
  end    = int(ipaddress.IPv4Address(iprange[1]))
  while(start < end):
    ips  = str(ipaddress.IPv4Address(start))
    start += 1
    if(not ping80(ips)):
      print "no response from:",  ips
      continue
    html = ""
    try:
      br   = initBrowser()
      r    = br.open("http://" + ips)
      html = r.read()
    except:
      print "mechanize error @ " + ips
      continue
    # saves.
    try:
      dir = outbase
      os.system("mkdir -p " + dir)
    except:
      print "mkdir error"
    try:
      f = open(dir + ips + "-inner.txt", "w")
      f.write(innerHTMLs(html))
      f.close()
      f = open(dir + ips + "-links.txt", "w")
      f.write(cutLinks(html))
      f.clos()
    except:
      print "file open error."

