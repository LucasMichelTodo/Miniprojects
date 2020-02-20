from time import time
import numpy as np
import multiprocessing as mp

print("Number of processors: ", mp.cpu_count())


# Prepare data: list of 200000 lists with 5 int each

np.random.RandomState(100)
arr = np.random.randint(0, 10, size=[200000, 5])
data = arr.tolist()
data[:5]

# Solution Without Paralleization


def howmany_within_range(row, minimum, maximum):
    """Returns how many numbers lie within `maximum` and `minimum` in a given `row`"""
    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
    return count


results = []
for row in data:
    results.append(howmany_within_range(row, minimum=4, maximum=8))

print(results[:10])


# Parallelizing using Pool.apply()

# Step 1: Init multiprocessing.Pool()
pool = mp.Pool(4)

# Step 2: `pool.apply` the `howmany_within_range()`
results = [pool.apply(howmany_within_range, args=(row, 4, 8)) for row in data]

# Step 3: Don't forget to close
pool.close()

print(results[:10])
# > [3, 1, 4, 4, 4, 2, 1, 1, 3, 3]
