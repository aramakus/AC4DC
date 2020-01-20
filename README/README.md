# AC4DC
Atomic/plasma code for XFEL damage simulations.

The code simulates the changes in electronic structure of atoms as they are illuminated by an XFEL pulse. Independent atom approximation enables to represent a molecule as collection of independent atoms, submerged into a bath of electrons, that are formed as a result of ionization. Main approximations:

1) All atoms interact with X-rays independent of each other.
2) At any time, an atom is described by an average-over-configuration state. 
3) Time evolution of atoms is described by a system of coupled rate equation (perturbation theory).
4) Electron plasma is split into two components - energetic photo-electrons and cold secondary electrons.
5) Photo-electrons have delta-function distribution around time-dependent average kinetic energy. No three-body reconbination (TBR) for photo-electrons, which can escape the sample.
6) Secondary electrons have Boltzmann distribution with time-dependent number density and temperature. They are assumed to be trapped at all times. Electron impact ionization (EII) and TBR are determined by average number density of secondary electrons over the volume of a sample.
