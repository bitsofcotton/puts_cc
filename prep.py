#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-

import sys

reload(sys)
sys.setdefaultencoding("utf-8")

import _mysql

conn = _mysql.connect(user='root', passwd='', host='localhost', db='wiktionary')
conn.query("select page_latest from page where page_title='" + _mysql.escape_string(sys.argv[1]) + "'")
pages = conn.use_result()
res   = ""
for s in range(0, 20):
  try:
    conn2 = _mysql.connect(user='root', passwd='', host='localhost', db='wiktionary')
    conn2.query("select old_text from text where old_id='" + _mysql.escape_string( pages.fetch_row()[0][0]) + "'")
    res += conn2.use_result().fetch_row()[0][0]
  except:
    break
print res

