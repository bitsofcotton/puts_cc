import sys

if(sys.argv[1] == '-'):
  sys.stdin.readline()
  sys.stdin.readline()
  b = int(sys.stdin.readline())
  s   = []
  B   = 0
  tb  = 1
  ctr = 0
  for line in sys.stdin:
    if(ctr == 0):
      B  = 0
      tb = 1
    ctr += 1
    if(int(line) >= b / 2): B += tb
    tb *= 2
    if(tb >= 255): ctr = 0
    if(ctr == 0):
      s.append(int(B))
  print(bytearray(s).decode("utf-8", errors="ignore"))
elif(sys.argv[1] == '+'):
  s = ""
  for line in sys.stdin:
    s += line[:- 1]
  s = s.encode("utf-8")
  print("P2")
  print("8", len(s))
  print("255")
  for ss in s:
    tb = 1
    for tt in range(0, 8):
      print((int(ss / tb) % 2) * 255)
      tb *= 2

