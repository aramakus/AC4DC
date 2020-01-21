# AC4DC
Atomic/plasma code for XFEL damage simulations.

The code simulates the changes in electronic structure of a atoms in bulk of a micro-crystal or amorphous media as it is illuminated by an X-ray free-electron laser (XFEL) pulse. An independent atom approximation is adopted to represent a molecule as a collection of independent atoms, submerged into a bath of electrons that are formed as a result of ionization.

### Installation

Compile script "compile_script.sh" is located in the script folder. Change 'g++' (linux, Ubuntu) to your compiler if required. The line "#include <eigen3/Eigen/Dense>" in the header file "EigenSolver.h" should reflect the location of the Eigen 3 library.  

### Dependencies:

+ Eigen 3        (vector algebra library, http://eigen.tuxfamily.org/index.php?title=Main_Page)
+ Wigner Symbols (3-j and 6-j symbols, https://github.com/joeydumont/wignerSymbols)
+ OpenMP         (multi-threading, https://www.openmp.org/)

### Key Approximations:

1) All atoms interact with X-rays independent of each other.
2) At any time, an atom is described by an average-over-configuration state. 
3) Time evolution of atoms is described by a system of coupled rate equation (perturbation theory).
4) Electron plasma is split into two components - energetic photo-electrons and cold secondary electrons.
5) Photo-electrons have delta-function distribution around time-dependent average kinetic energy. No three-body reconbination (TBR) for photo-electrons, which can escape the sample.
6) Secondary electrons have Boltzmann distribution with time-dependent number density and temperature. They are assumed to be trapped at all times. Electron impact ionization (EII) and TBR are determined by average number density of secondary electrons over the volume of a sample.

### Bibliography:

+ A. Kozlov and H. M. Quiney, _Comparison of Hartree-Fock and local density exchange approximations for calculation of radiation damage dynamics of light and heavy atoms in the field of x-ray free electron laser_, Phys. Scripta **94**, 075404 (2019). DOI: [10.1088/1402-4896/ab097c](https://doi.org/10.1088/1402-4896/ab097c)
+ C. P. Bhalla, N. O. Folland, and M. A. Hein, _Theoretical K-Shell Auger Rates, Transition Energies, and Fluorescence Yields for Multiply Ionized Neon_, Phys. Rev. A **8**, 649 (1973). DOI: [10.1103/PhysRevA.8.649](https://doi.org/10.1103/PhysRevA.8.649)
+ O. Yu. Gorobtsov, U. Lorenz, N. M. Kabachnik, and I. A. Vartanyants, _Theoretical study of electronic damage in single-particle imaging experiments at x-ray free-electron lasers for pulse durations from 0.1 to 10 fs_, Phys. Rev. E **91**, 062712 (2015). DOI: [10.1103/PhysRevE.91.062712](https://doi.org/10.1103/PhysRevE.91.062712)
+ W. R. Johnson, _Atomic Structure Theory: Lectures on Atomic Physics_, (Springer, Berlin, Heidelberg, 2007). DOI: [10.1007/978-3-540-68013-0](https://doi.org/10.1007/978-3-540-68013-0)