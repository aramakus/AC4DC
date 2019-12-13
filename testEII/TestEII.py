import matplotlib.pyplot as plt
import numpy as np
from math import log

File = open("EII.txt")
# read the content into a list "Rates"
data = []
for line in File:
    data.append(line.split())
    for i, elem in enumerate(data[-1]):
        data[-1][i] = float(elem)
File.close()
data = np.array(data).transpose()


plt.grid(True, which="both")
plt.title("Electron impact ionization of C(2p^-1)")
plt.xlabel("Electron energy, eV")
plt.ylabel("Cross-section, Mbarns")
plt.semilogx(data[0], data[1], lw = 2, c = 'red', alpha = 0.7)

plt.show()