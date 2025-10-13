import sys

if(sys.argv[1] == '-'):
  sys.stdin.readline()
  sys.stdin.readline()
  b = int(sys.stdin.readline())
  s = []
  B = 0
  c = 0
  for line in sys.stdin:
    c += 1
    if(int(line) >= b / 2): B += 1
    if(c >= 8):
      s.append(int(B))
      B = c = 0
    B *= 2
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
    ls = ss
    for tt in range(0, 8):
      print(ls)
      ls *= 2
      if(ls >= 256):
        ls -= 256
        ls += 1

