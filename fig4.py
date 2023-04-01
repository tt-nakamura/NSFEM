import numpy as np
import matplotlib.pyplot as plt
from NSFEM import NSFEM
from contour import contour

ns = NSFEM('cylinder/node.txt',
           'cylinder/element.txt',
           'cylinder/stream_fix.txt',
           'cylinder/vortex_fix.txt',
           'cylinder/vortex_wall.txt')
region = np.loadtxt('cylinder/region.txt')

# Karman vortex appears behind the cylinder
n = ns.run(0.5, 0.001, 0.003, maxit=440)
print(n)

plt.figure(figsize=(5,8.2))

x,y = ns.node.T
X,Y = region.T

plt.subplot(311)
plt.axis('equal')
plt.axis([0,20,0,10])
ns.plot_velocity()
plt.plot([0,20,20,0,0], [0,0,10,10,0], 'k')
plt.plot(X, Y, 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

plt.subplot(312)
plt.axis('equal')
plt.axis([0,20,0,10])
contour(x, y, ns.stream, 24, exclude=region, cmap='jet')
plt.plot([0,20,20,0,0], [0,0,10,10,0], 'k')
plt.plot(X, Y, 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

plt.subplot(313)
plt.axis('equal')
plt.axis([0,20,0,10])
contour(x, y, ns.vortex, 32, exclude=region, cmap='jet')
plt.plot([0,20,20,0,0], [0,0,10,10,0], 'k')
plt.plot(X, Y, 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

plt.tight_layout()
plt.savefig('fig4.eps')
plt.show()
