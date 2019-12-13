#pragma once

#include <vector>
#include "CoupledFunction.h"

class RadialWF : public CoupledFunction
{
public:

	RadialWF(int size = 0) : CoupledFunction(size)
	{
		infinity  = 0;
		turn = 0;
		Energy = 0;
		occup_number = 0;
	}
	RadialWF(const RadialWF& other);

	double Energy = 0;

	int L() { return l; }
	void set_L(int X);

	int GetNodes() { return (n-l-1); }
	int check_nodes(); //checks current number of nodes in F. Used in many routines to adjust the energy.
	int N() { return n; }

	void set_N(int X) { n = X; }
	void set_infinity(int X) { infinity = X; }
	void set_turn(int X) { turn = X; }// calculate number of nodes till turning point, discard the nodes after the turning point
	void set_occupancy(int X) { occup_number = X; }//set custom occupancy, useful for average over configurations

	int pract_infinity() { return infinity; }
	int occupancy() { return occup_number; }//might be fractional for future uses
	int turn_pt() { return turn; }

	const RadialWF& operator=(const RadialWF& Psi)
	{
		CoupledFunction::operator=(Psi);

		l = Psi.l;
		n = Psi.n;
		Energy = Psi.Energy;
		infinity = Psi.infinity;
		turn = Psi.turn;
		occup_number = Psi.occup_number;

		return *this;
	}

	~RadialWF() {}

protected:
	int l;
	int n;
	int infinity;
	int turn;
	int occup_number;
};
