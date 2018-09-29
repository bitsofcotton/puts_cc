#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-

import sys

reload(sys)
sys.setdefaultencoding("utf-8")

import MySQLdb

conn = MySQLdb.connect(user='root', passwd='', host='127.0.0.1', db='wiktionary', charset='utf8')
cur  = conn.cursor()
cur.execute("select page_latest from page where page_title=%s", [sys.argv[1]])
res  = ""
for s in range(0, 20):
  try:
    conn2 = MySQLdb.connect(user='root', passwd='', host='127.0.0.1', db='wiktionary', charset='utf8')
    cur2 = conn2.cursor()
    work = cur.fetchone()[0]
    cur2.execute("select old_text from text where old_id=%s", [work])
    res += cur2.fetchone()[0]
  except:
    break
if(0 < len(res.strip())):
  print res.strip(), 

