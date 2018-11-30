#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-

import mechanize
import re
import ipaddress
import socket
import os
import datetime
import glob
import hashlib
import sys

reload(sys)
sys.setdefaultencoding("utf-8")
mechanize._sockettimeout._GLOBAL_DEFAULT_TIMEOUT = 100
textmime = re.compile(r"(\.(txt|html|htm|shtml|shtm|xml|mml)$|/[\/\.]*$)")
elimiter = r"[ \t]+"

ipranges = [[u"153.126.0.0", u"153.126.127.255"],
            [u"153.126.128.0", u"153.126.255.255"],
            [u"153.127.0.0", u"153.127.127.255"],
            [u"153.127.128.0", u"153.127.191.255"]
            ]

inbase   = '~/Sites/datas/*/urls.txt'
outbase  = '~/Downloads/'
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
  br.set_handle_robots(True)
  return br

def cutLinks(html):
  urlf  = urls.findall(html)
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
    return False
  if result != 0:
    return False
  return True

def writeSub(stra, f, flog):
  if(not textmime.search(stra)):
    return
  html = ""
  try:
    br   = initBrowser()
    r    = br.open(stra)
    html = r.read()
    flog.write(datetime.datetime.now().strftime("%c") + stra + " : " + hashlib.sha256(html).hexdigest() + "\n")
    f.write(re.sub(elimiter, '', innerHTMLs(html).decode("utf-8")))
  except Exception as inst:
    sinst = str(repr(inst))
    if(re.search(r"robots\.txt", sinst)):
       return
    print "mechanize err @ " + stra
  except:
    print "mechanize error @ " + stra
  return

inurls = glob.glob(inbase)
for addr in inurls:
  lbase = os.path.dirname(addr) + "/web/sub/" + datetime.datetime.now().strftime("%Y%m%d%H") + ".txt"
  rbase = os.path.dirname(addr) + "/webhash.txt"
  print lbase
  f     = open(addr, "r")
  strs  = f.readlines()
  f.close()
  flog  = open(rbase, "a")
  f     = open(lbase, "w")
  for str0 in strs:
    stra, num = str0.split(" ")
    if(int(num) > 0):
      html = ""
      try:
        br = initBrowser()
        r  = br.open(stra)
        li = []
        for l in br.links():
          li.append(l.url)
        li = list(set(li))
      except:
        print "mechanize error @ " + stra
        continue
      for strs1 in li:
        writeSub(strs1, f, flog)
    else:
      writeSub(stra, f, flog)
  flog.close()
  f.close()

lbase = outbase + "/" + datetime.datetime.now().strftime("%Y%m%d%H")
fn    = open(outbase + "count.txt", "r")
num   = fn.readlines()
fn.close()
cnt   = int(num[0])
num0  = cnt
num   = 100
numl  = num
f     = open(lbase + ".txt",     "w")
fl    = open(lbase + "link.txt", "w")
for iprange in ipranges:
  if(num <= 0):
    break
  start  = int(ipaddress.IPv4Address(iprange[0]))
  end    = int(ipaddress.IPv4Address(iprange[1]))
  if(end - start <= cnt):
    cnt -= end - start
    continue
  start += cnt
  cnt    = 0
  while(start < end and cnt < num):
    ips    = str(ipaddress.IPv4Address(start))
    print ips, cnt, num
    num   -= 1
    start += 1
    if(not ping80(ips)):
      f.write("\n\n" + ips + "\nno response from : " + ips)
      continue
    html = ""
    try:
      br   = initBrowser()
      r    = br.open("http://" + ips)
      html = r.read()
    except:
      print "mechanize error @ " + ips
      continue
    try:
      f.write("\n\n" + ips + "\n")
      f.write(re.sub(elimiter, '', innerHTMLs(html).decode("utf-8")))
      fl.write("\n\n" + ips + "\n")
      fl.write(cutLinks(html))
    except:
      print "file open error."
f.close()
fl.close()
fn = open(outbase + "count.txt", "w")
if(cnt < num):
  fn.write("0")
else:
  fn.write(str(num0 + numl))
fn.close()

