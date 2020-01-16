#pragma once

#include <vector>

class PairFunction
{
public:
	PairFunction(int size = 0)
	{
		clear();
		if (size > 0) {
			resize(size);
		}
	}

	std::vector<double> F, G;  // Large and small components of radial wavefunction and it's derivatives

	void clear()
	{
		F.clear();
		G.clear();
	}

	void resize(int size)
	{
		F.resize(size);
		G.resize(size);
	}

	int size() { return static_cast<int>(F.size()); }

	void scale(double norm)
	{
		for (int i = 0; i < size(); i++) {
			F[i] *= norm;
			G[i] *= norm;
		}
	}

	const PairFunction& operator=(const PairFunction& Psi)
	{
		F = Psi.F;
		G = Psi.G;

		return *this;
	}

	~PairFunction() {}
};
