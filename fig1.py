import numpy as np
import matplotlib.pyplot as plt 
from NSFEM import NSFEM

ns = NSFEM('obstacle/node.txt',
           'obstacle/element.txt')
X,Y = np.loadtxt('obstacle/region.txt').T

plt.figure(figsize=(5,3))

plt.subplot(211)
plt.axis('equal')
plt.axis([0,20,0,5])
plt.plot(X,Y,'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

plt.subplot(212)
plt.axis('equal')
plt.axis([0,20,0,5])
ns.plot_element('k', lw=1)
plt.text(0, 5, 'A', ha='right')
plt.text(0, 0, 'B', ha='right', va='top')
plt.text(3.33, 0, 'C', va='top')
plt.text(6.67, 0, 'D', ha='right', va='top')
plt.text(20, 0, 'E', va='top')
plt.text(20, 5, 'F')
plt.box('off')
plt.xticks(())
plt.yticks(())
plt.tight_layout()
plt.savefig('fig1.eps')
plt.show()
