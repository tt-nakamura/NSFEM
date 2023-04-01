import numpy as np
import matplotlib.pyplot as plt
from contour import contour

plt.figure(figsize=(5,5))

x,y,u,v = np.loadtxt('fig6uv.txt').T
plt.subplot(221)
plt.axis('equal')
plt.axis([0,10,0,10])
plt.quiver(x,y,u,v);
plt.plot([0,10,10,0,0], [0,0,10,10,0], 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

x,y,ps = np.loadtxt('fig6ps.txt').T
plt.subplot(222)
plt.axis('equal')
plt.axis([0,10,0,10])
contour(x, y, ps, 20)
plt.plot([0,10,10,0,0], [0,0,10,10,0], 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

x,y,om = np.loadtxt('fig6om.txt').T
plt.subplot(223)
plt.axis('equal')
plt.axis([0,10,0,10])
contour(x, y, om, 20)
plt.plot([0,10,10,0,0], [0,0,10,10,0], 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

plt.tight_layout()
plt.savefig('fig6.eps')
plt.show()
