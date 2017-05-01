import sys
import numpy
import random

# receive as input the number of numbers followed by the name of the distribution
# uniform distribution is 1 to 10 * N
# print sorted seq of numbers NL separated
# <N>
# <num1>
# ...
# <numN>

if len(sys.argv) < 3:
    print("Expecting gendata.py nNums numDistribution")
    sys.exit(0)

n = int(sys.argv[1])
dist = sys.argv[2]
if dist == "uniform":
    nums = [random.randint(1,2**63 - 1) for _ in range(n)]
elif dist == "zipf05" or dist == "zipf25":
    nums = numpy.random.zipf(1.05 if dist == "zipf05" else 1.25, n)
elif dist == "uniform7":
    nums = [random.randint(1, 7 * 2**60) for _ in range(n)]
elif dist == "uniform_even":
    nums = [random.randint(1, 0x5555555555555555) for _ in range(n)]
else:
    print("Unexpected distribution." + dist + " Expecting uniform|zipf05|zipf25")
    assert(False)
print(n)
print('\n'.join([str(x) for x in sorted(nums)]))
