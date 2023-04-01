import numpy as np
import matplotlib.pyplot as plt 
from NSFEM import NSFEM

ns = NSFEM('cylinder/node.txt',
           'cylinder/element.txt')
X,Y = np.loadtxt('cylinder/region.txt').T

plt.figure(figsize=(5,5.4))

plt.subplot(211)
plt.axis('equal')
plt.axis([0,20,0,10])
plt.plot(X,Y,'k')
plt.plot([0,20,20,0,0], [0,0,10,10,0], 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

plt.subplot(212)
plt.axis('equal')
plt.axis([0,20,0,10])
ns.plot_element('k', lw=1)
plt.text(0, 10, 'A', ha='right')
plt.text(0, 0, 'B', ha='right', va='top')
plt.text(20, 0, 'C', va='top')
plt.text(20, 10, 'D')
plt.box('off')
plt.xticks(())
plt.yticks(())
plt.tight_layout()
plt.savefig('fig3.eps')
plt.show()
