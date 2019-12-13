#include <vector>
#include <iostream>
#include <algorithm>
#include "Spline.h"

void spln_iter(std::vector<double>* , unsigned int, std::vector<double>*, double*, unsigned int);
void dspln_iter(std::vector<double>*, unsigned int, std::vector<double>*, double*, unsigned int);

Bspline::Bspline(std::vector<double>* Knots)
{
	knots = Knots;
	size = (*Knots).size() - 1;
}

Bspline::Bspline(std::vector<double>* Knots, int Order)
{
	knots = Knots;
	size = (*Knots).size() - 1;
	order = Order;
}


double Bspline::Value(unsigned int k, unsigned int i, double x)
{   
	double Result, A, B;
	order = k;

	if ( (x <= (*knots)[i]) || (x > (*knots)[i + k + 1]) || (i+k+1)>(*knots).size() || (i < 0) ) { Result = 0.; }
	else
	{
		if (k == 0)
		{
			Result = 1.;
		}
		else
		{
			if ((*knots)[i] == (*knots)[i + k]) { A = 0.; }
			else { A = (x - (*knots)[i]) / ((*knots)[i + k] - (*knots)[i]); }
			if ((*knots)[i + 1] == (*knots)[i + k + 1]) { B = 0.; }
			else { B = ((*knots)[i + k + 1] - x) / ((*knots)[i + k + 1] - (*knots)[i + 1]); }
			Result = Value(k - 1, i, x)*A + Value(k - 1, i + 1, x)*B;
		}
	}

	return Result;
}

double Bspline::dValue(unsigned int k, unsigned int i, unsigned int d_Order, double x)
{
	double Result, A, B;
	order = k;

	if ((x <= (*knots)[i]) || (x > (*knots)[i + k + 1]) || (i + k + 1)>(*knots).size() || (i < 0)) { Result = 0.; }
	else
	{
		if (d_Order > k)
		{
			Result = 0.;
		}
		else
		{
			if (d_Order == 0)
			{
				Result = Value(k, i, x);
			}
			else
			{
				if ((*knots)[i] == (*knots)[i + k]) { A = 0.; }
				else { A = (double)(k) / ((*knots)[i + k] - (*knots)[i]); }
				if ((*knots)[i + 1] == (*knots)[i + k + 1]) { B = 0.; }
				else { B = (double)(k) / ((*knots)[i + k + 1] - (*knots)[i + 1]); }
				Result = dValue(k - 1, i, d_Order - 1, x)*A - dValue(k - 1, i + 1, d_Order - 1, x)*B;
			}
		}
	}


	return Result;
}

std::vector<double> Bspline::Values(unsigned int k, unsigned int d_Order, double x)
{
	order = k;

	unsigned int Right = k + 1; // closest to "x" knot from the right hand side.
	int NumSplines = 0; // number of non-trivial splines at "x"
	unsigned int tmp_k = 0, tmp_d_order = 0;
	values.clear();
	spln_nums.clear();

	while (((*knots)[Right] < x) && (Right < ((*knots).size() - 1 - k)) )
	{
		Right++;
	}
	if ((*knots)[Right] < x) { std::cout << "attemt to calculate Bspline value outside of knots range \n";  }
	else
	{
		//find number of non-trivial splines at "x"
		NumSplines = std::min(k, ((*knots).size() - Right - 1)) + std::min(k, (Right - 1)) - k + 1;
		if (NumSplines < 0) { NumSplines = 0; }

		else
		{
			values.resize(NumSplines);
			values[NumSplines - 1] = 1.;
			spln_nums.resize(NumSplines);
			values[NumSplines-1] = 1;
			for (unsigned int i = NumSplines; i > 0; i--)
			{
				spln_nums[i-1] = Right - NumSplines  + i - 1;
			}
			
			while (tmp_k < (k - d_Order))
			{
				spln_iter(&values, Right, knots, &x, tmp_k);
				tmp_k++;
			}

			while (d_Order > tmp_d_order)
			{
				dspln_iter(&values, Right, knots, &x, tmp_k);
				tmp_k++;
				tmp_d_order++;
			}
		}

	}

	return values;
}

double Bspline::Knot(int i)
{
	if (i < (*knots).size())
	{
		return (*knots)[i];
	}
	else { return -1.; }
}




void spln_iter(std::vector<double>* values, unsigned int right_knt, std::vector<double>* knots, double* x, unsigned int curr_k)
{
	// this function is used for iterative calculation of all non-trivial B-splines at "x"
	// it takes array of values, where previous iteration results are (B_{curr_k}),
	// which is stored in "values". If curr_k < k, first few elements of "values" are 0.
	// Routine calculates (B_{curr_k+1}) and puts back into *values. Number of first null elements in "values" is reduced by 1 as a result.
	double A, B;
	int TotNumSpl = (*values).size();
	int i = TotNumSpl - 1;

	std::vector<double> Tmp;
	Tmp.resize(TotNumSpl);


	A = (*x - (*knots)[right_knt - 1]) / ((*knots)[right_knt + curr_k] - (*knots)[right_knt - 1]);
	Tmp[i] = A*(*values)[i];
	i--;

	while ( (i >= 0) )
	{
		//somewhere here mistake. Check
		if ((*knots)[right_knt - TotNumSpl + i] == (*knots)[right_knt - TotNumSpl + i + curr_k + 1]) { A = 0.; }
		else { A = (*x - (*knots)[right_knt - TotNumSpl + i]) / ( (*knots)[right_knt - TotNumSpl + i + curr_k + 1] - (*knots)[right_knt - TotNumSpl + i] ); }
		if ((*knots)[right_knt - TotNumSpl + i + 1] == (*knots)[right_knt - TotNumSpl + i + curr_k + 2]) { B = 0.; }
		else { B = ((*knots)[right_knt - TotNumSpl + i + curr_k + 2] - *x) / ((*knots)[right_knt - TotNumSpl + i + curr_k + 2] - (*knots)[right_knt - TotNumSpl + i + 1]); }
		Tmp[i] = A*(*values)[i] + B*(*values)[i + 1];
		i--;
	}
	(*values) = Tmp;
}



void dspln_iter(std::vector<double>* values, unsigned int right_knt, std::vector<double>* knots, double* x, unsigned int curr_k)
{
	// similar to spln_iter,
	// this function is used for iterative calculation of all non-trivial B-splines Derivatives at "x"
	// it produces array of derivatives of (B_{curr_k+1}) out of (B_{curr_k}).
	// Routine calculates derivatives of (B_{curr_k+1}) and puts back into *values. Number of null elements in "values" is reduced by 1.
	double A, B;
	int TotNumSpl = (*values).size();
	int i = TotNumSpl - 1;

	std::vector<double> Tmp;
	Tmp.resize(TotNumSpl);


	A = (double)(curr_k + 1) / ((*knots)[right_knt + curr_k] - (*knots)[right_knt - 1]);
	Tmp[i] = A*(*values)[i];
	i--;

	while ((i >= 0) && ((*values)[i + 1] != 0))
	{
		//somewhere here mistake. Check
		if ((*knots)[right_knt - TotNumSpl + i] == (*knots)[right_knt - TotNumSpl + i + curr_k + 1]) { A = 0.; }
		else { A = (double)(curr_k + 1) / ((*knots)[right_knt - TotNumSpl + i + curr_k + 1] - (*knots)[right_knt - TotNumSpl + i]); }
		if ((*knots)[right_knt - TotNumSpl + i + 1] == (*knots)[right_knt - TotNumSpl + i + curr_k + 2]) { B = 0.; }
		else { B = (double)(curr_k + 1) / ((*knots)[right_knt - TotNumSpl + i + curr_k + 2] - (*knots)[right_knt - TotNumSpl + i + 1]); }
		Tmp[i] = A*(*values)[i] - B*(*values)[i + 1];
		i--;
	}
	(*values) = Tmp;
}