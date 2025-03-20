import os
import sys

start = int(sys.argv[1])
end = int(sys.argv[2])
suffix = sys.argv[3]

should_have = list(range(start, end + 1))

for file in os.listdir():
    if suffix in file:
        should_have.remove(int(file.replace("Run", "").replace(suffix, "")))
print(should_have)
