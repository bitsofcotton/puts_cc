#! /usr/bin/env python

import os
import sys
import glob
import subprocess

dbase   = '~/Downloads/texts/crawled'
prog    = '~/Sites/puts'
inbase  = '~/Sites/datas/cd6d67d2d39b15de942baa6df6f6144522042d2d7e7acd9de4826e0809312665'
outbase = '~/Downloads/outputs'

dicts   = glob.glob(inbase + '/dicts/*')
topics  = glob.glob(inbase + '/topics/*')
words   = inbase + '/words.txt'

for r in glob.glob(dbase + '/*-inner.txt'):
  cmdl = 'cd ' + dbase + '/../output && ' + 'cat ' + r + ' | ' + prog + ' lword > ' + os.path.basename(r) + '-lword.txt'
  cmdc = 'cd ' + dbase + '/../output && ' + 'cat ' + r + ' | ' + prog + ' toc ' + words + ' ' + ' '.join(dicts) + ' -toc ' + ' '.join(topics) + ' > ' + os.path.basename(r) + '-corpus.txt'
  print cmdl, cmdc
  if(not os.path.exists(os.path.basename(r) + '-lword.txt')):
    subprocess.call(cmdl, shell=True)
  if(not os.path.exists(os.path.basename(r) + '-corpus.txt')):
    subprocess.call(cmdc, shell=True)

