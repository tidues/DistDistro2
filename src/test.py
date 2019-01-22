import time

cnt = 16450

v1 = 0
v2 = 0
v3 = 0

start = time.time()
for i in range(cnt):
    for j in range(cnt):
        v1 += i + j
        v2 += i + j
        v3 += i + j
end = time.time()
print('nested: ', end - start)

start = time.time()
for i in range(cnt):
    for j in range(cnt):
        v1 += i + j

for i in range(cnt):
    for j in range(cnt):
        v2 += i + j

for i in range(cnt):
    for j in range(cnt):
        v3 += i + j
end = time.time()
print('splited: ', end - start)
        
