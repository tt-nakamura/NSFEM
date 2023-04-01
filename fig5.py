import numpy as np
import matplotlib.pyplot as plt 
from NSFEM import NSFEM

ns = NSFEM('cavity/node.txt',
           'cavity/element.txt')

plt.figure(figsize=(5,2.5))

plt.subplot(121)
plt.axis('equal')
plt.axis([0,10,0,10])
plt.plot([0,10,10,0,0], [0,0,10,10,0], 'k')
plt.box('off')
plt.xticks(())
plt.yticks(())

plt.subplot(122)
plt.axis('equal')
plt.axis([0,10,0,10])
ns.plot_element('k')
plt.text(0, 10, 'A', ha='right')
plt.text(0, 0, 'B', ha='right', va='top')
plt.text(10, 0, 'C', va='top')
plt.text(10, 10, 'D')
plt.box('off')
plt.xticks(())
plt.yticks(())
plt.tight_layout()
plt.savefig('fig5.eps')
plt.show()
