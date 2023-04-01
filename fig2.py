import numpy as np
import matplotlib.pyplot as plt
from NSFEM import NSFEM
from contour import contour

ns = NSFEM('obstacle/node.txt',
           'obstacle/element.txt',
           'obstacle/stream_fix.txt',
           'obstacle/vortex_fix.txt',
           'obstacle/vortex_wall.txt')
region = np.loadtxt('obstacle/region.txt')

n = ns.run(0.2, 0.05, 0.01)
print(n)

plt.figure(figsize=(5,4.5))

x,y = ns.node.T
X,Y = region.T

plt.subplot(311)
plt.axis('equal')
plt.axis([0,20,0,5])
ns.plot_velocity()
plt.plot(X, Y, 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

plt.subplot(312)
plt.axis('equal')
plt.axis([0,20,0,5])
contour(x, y, ns.stream, 16, confine=region, cmap='jet')
plt.plot(X, Y, 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

plt.subplot(313)
plt.axis('equal')
plt.axis([0,20,0,5])
contour(x, y, ns.vortex, 32, confine=region, cmap='jet')
plt.plot(X, Y, 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

plt.tight_layout()
plt.savefig('fig2.eps')
plt.show()
