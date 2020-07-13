#!/bin/zsh
#g++-9 installed via Homebrew. Standard Apple gcc distribution does not support openMP
g++-9 -O3 -std=c++11 -I"$HOME/Programming/include" main.cpp Constant.cpp DecayRates.cpp EigenSolver.cpp Grid.cpp HartreeFock.cpp Input.cpp IntegrateRateEquation.cpp Numerics.cpp Potential.cpp RadialWF.cpp RateEquationSolver.cpp Plasma.cpp Wigner/wignerSymbols.cpp -o HF -fopenmp
