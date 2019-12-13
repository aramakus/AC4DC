#include "SplineIntegrals.h"
#include <fstream>
#include <iomanip>
#include <math.h>

BsplineIntegrals::BsplineIntegrals(Bspline &Bas)
{
	double tmp;
	std::vector<double> B, d2_B;
	std::vector<int> NonZeroSplines;
	GaussLegendre I(8);

	overlap.clear();
	overlap.resize(Bas.Size());
	kinetic.clear();
	kinetic.resize(Bas.Size());
	hydro_potential.clear();
	hydro_potential.resize(Bas.Size());

	for (int i = 0; i < Bas.Size(); i++)
	{
		overlap[i].resize(Bas.Size());
		kinetic[i].resize(Bas.Size());
		hydro_potential[i].resize(Bas.Size());
	}

	for (int i = 0; i < (Bas.Size() + Bas.Order() - 1); i++)//total number of knots
	{
		if (Bas.Knot(i) != Bas.Knot(i + 1))
		{
			I.SetInterval(Bas.Knot(i), Bas.Knot(i + 1));
			for (int j = 0; j < I.GetOrder(); j++)
			{
				B = Bas.Values(Bas.Order(), 0, I.GetAbscissa(j));
				NonZeroSplines = Bas.Values_i();
				d2_B = Bas.Values(Bas.Order(), 2, I.GetAbscissa(j));
				NonZeroSplines = Bas.Values_i();

				for (int n = 0; n < NonZeroSplines.size(); n++)
				{
					for (int m = n; m < NonZeroSplines.size(); m++)
					{
							tmp = I.GetWeight(j) * B[n] * B[m];
							overlap[NonZeroSplines[n]][NonZeroSplines[m]] += tmp;
							kinetic[NonZeroSplines[n]][NonZeroSplines[m]] -= 0.5*I.GetWeight(j) * B[n] * d2_B[m];
							hydro_potential[NonZeroSplines[n]][NonZeroSplines[m]] -= tmp / I.GetAbscissa(j) * (1. - 1. / I.GetAbscissa(j));// s-wave case. For non-zero orbital momentum extra term needed.
					}
				}
			}
		}
	}

	for (int i = 0; i < Bas.Size(); i++)
	{
		for (int j = 0; j < i; j++)
		{
			overlap[i][j] = overlap[j][i];
			kinetic[i][j] = kinetic[j][i];
			hydro_potential[i][j] = hydro_potential[j][i];
		}
	}

	std::ofstream FL1("overlapDebug.txt");
	FL1.precision(4);

	for (int i = 0; i < Bas.Size(); i++)
	{
		for (int j = 0; j < Bas.Size(); j++)
		{
			FL1 << std::setw(10) << kinetic[i][j] + hydro_potential[i][j] << " ";
		}
		FL1 << "\n";
	}
	FL1.close();


}