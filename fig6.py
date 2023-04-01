import numpy as np
import matplotlib.pyplot as plt
from NSFEM import NSFEM
from contour import contour

ns = NSFEM('cavity/node.txt',
           'cavity/element.txt',
           'cavity/stream_fix.txt',
           'cavity/vortex_fix.txt',
           'cavity/vortex_wall.txt')

n = ns.run(0.1, 0.01, 0.01)
print(n)

plt.figure(figsize=(5,5))

x,y = ns.node.T

plt.subplot(221)
plt.axis('equal')
plt.axis([0,10,0,10])
ns.plot_velocity()
plt.plot([0,10,10,0,0], [0,0,10,10,0], 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

plt.subplot(222)
plt.axis('equal')
plt.axis([0,10,0,10])
contour(x, y, ns.stream, 20)
plt.plot([0,10,10,0,0], [0,0,10,10,0], 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

plt.subplot(223)
plt.axis('equal')
plt.axis([0,10,0,10])
contour(x, y, ns.vortex, 20)
plt.plot([0,10,10,0,0], [0,0,10,10,0], 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

plt.tight_layout()
plt.savefig('fig6.eps')
plt.show()
