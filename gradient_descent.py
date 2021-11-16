import plotly.graph_objects as go
from random import randint

# Paraboloid function
def paraboloid(x, y) -> float:
    return x ** 2 + y ** 2

# Test data generation (only really necessary for the plotting below)
xs_start = ys_start = -10
xs_stop = ys_stop = 11
xs_step = ys_step = 1

xs = [i for i in range(xs_start, xs_stop, xs_step)]
ys = [i for i in range(ys_start, ys_stop, ys_step)]
zs = []

for x in xs:
    temp_res = []
    for y in ys:
        result = paraboloid(x, y)
        temp_res.append(result)
    zs.append(temp_res)

# Plot paraboloid
fig = go.Figure(go.Surface(x=xs, y=ys, z=zs, colorscale='Viridis'))
fig.show()

# Gradient descend
def compute_gradient(vec):
    assert len(vec) == 2
    x = vec[0]
    y = vec[1]
    return [2 * x, 2 * y]

def compute_step(curr_pos, learning_rate):
    grad = compute_gradient(curr_pos)
    grad[0] *= -learning_rate
    grad[1] *= -learning_rate
    next_pos = [0, 0]
    next_pos[0] = curr_pos[0] + grad[0]
    next_pos[1] = curr_pos[1] + grad[1]
    return next_pos

# Starting position
start_pos = []

# Ensure that we don't start at a minimum (0, 0 in our case)
while True:
    start_x = randint(xs_start, xs_stop)
    start_y = randint(ys_start, ys_stop)
    if start_x != 0 and start_y != 0:
        start_pos = [start_x, start_y]
        break

print(start_pos)

# Run gradient descend
epochs: int = 5000
learning_rate = 0.001
    
best_pos = start_pos

for i in range(0, epochs):
    next_pos = compute_step(best_pos, learning_rate)
    if i % 500 == 0:
        print(f'Epoch {i}: {next_pos}')
    best_pos = next_pos    

print(f'Best guess for a minimum: {best_pos}')

