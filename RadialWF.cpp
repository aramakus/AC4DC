#include "RadialWF.h"

RadialWF::RadialWF(const RadialWF& other) : CoupledFunction(other), l(other.l), Energy(other.Energy), n(other.n), infinity(other.infinity), turn(other.turn), occup_number(other.occup_number)
{}

int RadialWF::check_nodes()
{
	int Result = 0;
	for (int i = 1; i < turn; i++) {
		if ((F[i] * F[i - 1]) < 0) {
			Result++;
		}
	}

	return Result;
}

void RadialWF::set_L(int X)
{
	l = X;
	occup_number = 4 * l + 2;
}
	/* file Read/Write should be added here later, following Julz.*/
