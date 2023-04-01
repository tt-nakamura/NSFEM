import numpy as np
import matplotlib.pyplot as plt
from contour import contour

region = np.loadtxt('../obstacle/region.txt')
X,Y = region.T

plt.figure(figsize=(5, 4.5))

x,y,u,v = np.loadtxt('fig2uv.txt').T
plt.subplot(311)
plt.axis('equal')
plt.axis([0,20,0,5])
plt.quiver(x,y,u,v);
plt.plot(X, Y, 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

x,y,ps = np.loadtxt('fig2ps.txt').T
plt.subplot(312)
plt.axis('equal')
plt.axis([0,20,0,5])
contour(x, y, ps, 16, confine=region, cmap='jet')
plt.plot(X, Y, 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

x,y,om = np.loadtxt('fig2om.txt').T
plt.subplot(313)
plt.axis('equal')
plt.axis([0,20,0,5])
contour(x, y, om, 16, confine=region, cmap='jet')
plt.plot(X, Y, 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

plt.tight_layout()
plt.savefig('fig2.eps')
plt.show()
