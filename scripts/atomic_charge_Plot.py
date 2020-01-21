import matplotlib.pyplot as plt
import numpy as np
from math import log
import os.path as path

path = path.abspath(path.join(__file__ ,"../../output"))

File = open(path + "/Charge_O.txt")
# read the content into a list "Rates"
charge = []
for line in File:
    charge.append(line.split())
    for i, elem in enumerate(charge[-1]):
        charge[-1][i] = float(elem)
File.close()
charge = charge[100:-100]
charge = np.array(charge)
charge = charge.transpose()

File = open(path + "/Intensity_water.txt")
# read the content into a list "Rates"
intensity = []
for line in File:
    intensity.append(line.split())
    for i, elem in enumerate(intensity[-1]):
        intensity[-1][i] = float(elem)
File.close()
intensity = intensity[100:-100]
intensity = np.array(intensity)
intensity = intensity.transpose()

plt.rcParams["font.size"] = 14
fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(111)

for i, elem in enumerate(charge[:-1]):
    ax.plot(charge[-1], elem, lw = 2, alpha = 0.7, label = str(i))

ax.plot(intensity[1], intensity[0], lw = 2, c = 'black', ls = '--', alpha = 0.7)
ax.set_title("Charge state dynamics")
ax.set_xlabel("Time, fs")
ax.set_ylabel("Probability")

plt.figlegend(loc = (0.11, 0.43))
plt.subplots_adjust(left=0.1, right=0.92, top=0.93, bottom=0.1)

plt.show()