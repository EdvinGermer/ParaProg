import numpy as np

with open('/home/maya/public/PDP_Assignment1/input96.txt', 'r') as f:
    fx = np.fromstring(f.read(), sep=' ')
    fx = fx[1:]


    fx = fx.reshape(8, 12)

    np.savetxt('test.txt', fx, delimiter=',', fmt='%f')