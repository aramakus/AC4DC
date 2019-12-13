#include "EigenSolver.h"
#include <vector>
#include "stdafx.h"

EigenSolver::EigenSolver()
{
	//this class would test sever different pakages for matrix diagonalization. It uses just "Eigen" for now. Reminder: "Armadillo" in future
}

void EigenSolver::SolveSystem(std::vector<std::vector<double>> &A, std::vector<double> &B, int Length)
{
	Eigen::MatrixXd M(Length, Length);
	Eigen::VectorXd Y(Length);
	Eigen::VectorXd X(Length);

	for (int i = 0; i < Length; i++)
	{
		Y(i) = B[i];
		for (int j = 0; j < Length; j++)
		{
			M(i, j) = A[i][j];
		}
	}

	X = M.fullPivLu().solve(Y);

	for (int i = 0; i < Length; i++)
	{
		B[i] = X(i);
	}
}

void EigenSolver::SolveGenEig(std::vector<std::vector<double>> &F, std::vector<std::vector<double>> &S)
{
	int Length = F.size();
	Eigen::MatrixXd M(Length, Length);
	Eigen::MatrixXd Q(Length, Length);

	for (int i = 0; i < Length; i++) {
		for (int j = 0; j < Length; j++) {
			M(i, j) = F[i][j];
			Q(i, j) = S[i][j];
		}
	}

	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> Magic(M, Q);

	EigenValues = Magic.eigenvalues();
	EigenVectors = Magic.eigenvectors();
}

vector<double> EigenSolver::EigenVals()
{
	vector<double> Result(EigenValues.size(), 0);

	for (int i = 0; i < Result.size(); i++) {
		Result[i] = EigenValues(i);
	}

	return move(Result);
}

vector<vector<double>> EigenSolver::EigenVecs()
{
	vector<vector<double>> Result(EigenValues.size(), vector<double>(EigenValues.size(), 0));
	for (int i = 0; i < Result.size(); i++) {
		for (int j = 0; j < Result.size(); j++) {
			Result[i][j] = EigenVectors(j, i);
		}
	}

	return move(Result);
}