import time
from sympy import *
from sympy.abc import x

cnt = 100000

begin = time.time()
t = 0
for i in range(cnt):
    res = t + 1
end = time.time()
print('normal:', end - begin)
    
begin = time.time()
for i in range(cnt):
    res = x + 1
end = time.time()
print('sympy:', end - begin)


