import matplotlib.pyplot as plt
import numpy as np
from math import log
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import os.path as path

chem_elems = {'1' : 'H', '6' : 'C', '7' : 'N', '8' : 'O', '15' : 'P', '16' : 'S'}

path = path.abspath(path.join(__file__ ,"../../output"))
File = open(path + "/MD_Data.txt")
# read the content into a list "Rates"
file_content = []
list_of_atoms = []
for ind, line in enumerate(File):
    if ind == 1:
        list_of_atoms = line.split()[1:-2]
        list_of_atoms = [chem_elems[n] for n in list_of_atoms]

    if ind < 2: continue
    file_content.append(line.split())
    for i, elem in enumerate(file_content[-1]):
        file_content[-1][i] = float(elem)
File.close()
file_content = np.array(file_content)

time = file_content[:, 0]
charge = file_content[:, 1:-2]
electrons = file_content[:, -2:]


# Color them according to the maximum charge.
max_chrg = charge[-1, :]

norm_max_chrg = (max_chrg / np.max( max_chrg ) + 0.1)/1.1

cNorm = colors.Normalize(vmin=0, vmax=1)
cm = plt.get_cmap('gist_rainbow')
scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)

File = open(path + "/Intensity_LassoPeptide.txt")
# read the content into a list "Rates"
intsty = []
for line in File:
  intsty.append(line.split())
  for i, elem in enumerate(intsty[-1]):
    intsty[-1][i] = float(elem)
File.close()
intsty = np.array(intsty)
intsty = intsty[:, 0]

electrons[0, :] = electrons[1, :]
Temperature = np.divide(electrons[:, 1], electrons[:, 0])*2/3 * 27.2

scaling = np.max(charge)*1.2

plt.rcParams["font.size"] = 14
fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(111)

for i, elem in enumerate(norm_max_chrg):
  ax.plot(time, charge[:, i], lw = 2, c = scalarMap.to_rgba(norm_max_chrg[i]), alpha = 0.7, label = list_of_atoms[i])

ax.plot(time, intsty*0.9*scaling, lw = 2, c = 'black', ls = '--', alpha = 0.7, label = "Pulse profile")
ax.set_xlabel('time (fs)')
ax.set_title('Secondary electrons temperature and ion charges')
ax.set_ylabel('Ionization degree')
plt.ylim(bottom=0, top = scaling)

ax2 = ax.twinx()
ax2.plot(time, Temperature, lw = 1, c = 'red', alpha = 0.7, label = "Plasma Temperature")
ax2.set_ylim(bottom = 0, top = Temperature[-1]*1.4)
ax2.set_ylabel('Temperature, eV')

plt.figlegend(loc = (0.11, 0.56))
plt.subplots_adjust(left=0.1, right=0.92, top=0.93, bottom=0.1)

plt.show()