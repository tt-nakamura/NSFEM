import numpy as np
import matplotlib.pyplot as plt
from contour import contour

region = np.loadtxt('../cylinder/region.txt')
X,Y = region.T

plt.figure(figsize=(5, 8.2))

x,y,u,v = np.loadtxt('fig4uv.txt').T
plt.subplot(311)
plt.axis('equal')
plt.axis([0,20,0,10])
plt.quiver(x,y,u,v);
plt.plot([0,20,20,0,0], [0,0,10,10,0], 'k')
plt.plot(X, Y, 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

x,y,ps = np.loadtxt('fig4ps.txt').T
plt.subplot(312)
plt.axis('equal')
plt.axis([0,20,0,10])
contour(x, y, ps, 20, exclude=region, cmap='jet')
plt.plot([0,20,20,0,0], [0,0,10,10,0], 'k')
plt.plot(X, Y, 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

x,y,om = np.loadtxt('fig4om.txt').T
plt.subplot(313)
plt.axis('equal')
plt.axis([0,20,0,10])
contour(x, y, om, 20, exclude=region, cmap='jet')
plt.plot([0,20,20,0,0], [0,0,10,10,0], 'k')
plt.plot(X, Y, 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

plt.tight_layout()
plt.savefig('fig4.eps')
plt.show()
