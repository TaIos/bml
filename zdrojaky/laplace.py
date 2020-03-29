import numpy as np
from scipy.stats import norm
import matplotlib.pylab as plt

def pdf(x):
    out = 0.8 * norm.pdf(x, loc=4, scale=4) 
    out += 0.1 * norm.pdf(x, loc=7,scale=1.3)
    out += 0.1 * norm.pdf(x, loc=11,scale=1.3)
    return out


x = np.linspace(-10, 18, 1000)

plt.figure(1, figsize=(14,4))
plt.subplot(1,3,1)
plt.plot(x, pdf(x))
plt.subplot(1,3,2)
plt.plot(x, pdf(x))
plt.axvline(x=6.4, color='black', lw=2, ls='--')
#plt.plot(x, norm.pdf(x, loc=6.5, scale=3), lw=3)
plt.subplot(1,3,3)
plt.axvline(x=6.4, color='black', lw=2, ls='--')
plt.plot(x, pdf(x))
plt.plot(x, norm.pdf(x, loc=6.4, scale=4.2), lw=3)
plt.savefig('/tmp/laplace.png', bbox_inches='tight')