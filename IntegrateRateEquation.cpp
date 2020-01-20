/*===========================================================================
This file is part of AC4DC.

    AC4DC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AC4DC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AC4DC.  If not, see <https://www.gnu.org/licenses/>.
===========================================================================*/
#include "IntegrateRateEquation.h"
#include <algorithm>
#include <iostream>
#include <omp.h>
#include <cmath>

static const double Moulton_5[5] = { 251. / 720., 646. / 720., -264. / 720., 106. / 720., -19. / 720. }; //Adams�Moulton method
static const double Bashforth_5[5] = { 1901. / 720., -1378. / 360., 109. / 30., -637. / 360., 251. / 720. }; //Adams�Bashforth method

inline bool CompareChar(vector<char>&, char);


IntegrateRateEquation::IntegrateRateEquation(vector<double> &dT, vector<double> &T, AtomRateData& Store, vector<double> InitCond, const vector<double>& Intensity) : dt(dT), t(T), f(Intensity), store(Store)
{
	adams_n = 5;
	t.resize(dt.size());
	if (f.size() != dt.size()) f = vector<double>(dt.size(), 0);
	dpdt.resize(InitCond.size());
	p.resize(InitCond.size());
	for (int i = 0; i < InitCond.size(); i++)
	{
		dpdt[i].resize(adams_n + 1);
		p[i].resize(adams_n + 1);
		p[i][0] = InitCond[i];
	}

	// Initialize p, dpdt for the subsequent integration
	vector<double> A(p.size(), 0);
	vector<double> X(p.size(), 0);
	double tollerance = 0.000001, error = 1, old_p = 0, tmp = 0;
	for (int m = 0; m < adams_n; m++) {
		if (m > 0) {
			for (auto& v : p) {// Guess.
				v[m] = v[m - 1];
			}
		}
		while (error > tollerance) {
			error = 0;

			for (auto& rate: Store.Photo) {
				tmp = rate.val*f[m];
				A[rate.from] -= tmp;
				X[rate.to] += tmp * p[rate.from][m];
			}
			for (auto& rate: Store.Fluor) {
				A[rate.from] -= rate.val;
				X[rate.to] += rate.val * p[rate.from][m];
			}
			for (auto& rate: Store.Auger) {
				A[rate.from] -= rate.val;
				X[rate.to] += rate.val * p[rate.from][m];
			}
			/*
			for (int i = 0; i < Store.Rates.size(); i++) {
				if (CompareChar(depends_t, Store.Rates[i].type)) {
					A[Store.Rates[i].from] -= Store.Rates[i].val*f[m];
					X[Store.Rates[i].to] += Store.Rates[i].val*f[m] * p[Store.Rates[i].from][m];
				} else {
					A[Store.Rates[i].from] -= Store.Rates[i].val;
					X[Store.Rates[i].to] += Store.Rates[i].val * p[Store.Rates[i].from][m];
				}
			}
			*/
			for (int i = 0; i < p.size(); i++)// Correct.
			{
				if (A[i] == 0 && X[i] == 0) continue;
				dpdt[i][m] = A[i]*p[i][m] + X[i];
				A[i] = 0.;
				X[i] = 0.;
				if (m > 0)
				{
					old_p = p[i][m];
					p[i][m] = p[i][m - 1] + dt[m - 1] * 0.5*(dpdt[i][m] + dpdt[i][m-1]);
					tmp = fabs(p[i][m] - old_p) / p[i][m];
					if (error < tmp) error = tmp;
				}
			}
		}
		error = 1;
	}
}


IntegrateRateEquation::IntegrateRateEquation(vector<double> &dT, vector<double> &T, vector<AtomRateData> & Store, Plasma & Elecs, const vector<double>& Intensity) : dt(dT), t(T), f(Intensity), store(Store[0])
{
	// f(F) is intensity defined at times T[m].
	// InitCond defines number of states and initial values for p.

	// Recast plasma equations (instant thermalization of photoelectrons)
	//
	//  dNp[t]/dt = n*Sum{i} p[i][t]*( Pht[i] - Gesc[t]*Np[t] )
	//  dEp[t]/dt = n*Sum{i} p[i][t]*( ePht[i] - Np[t]*(W[i] + eSp[i]) - Gesc[t]*Ep[t] - Cll[t](N,Np,E,Ep) )	
	//  dN[t]/dt = n*Sum{i} p[i][t]*( Aug[i][t] + N[t]*S[i][t] + Np[t]*Sp[i][t])
	//  dE[t]/dt = n*Sum{i} p[i][t]*( eAug[i][t] - N[t]*eS[i] + Np[t]*W[i])	
	//
	//	e[ij] - ionization potential of electron with ion going from "i" to "j"
	//  S[i][t] = Sum{j} * Sum{j} Int{_0^infty} dp*p^3*f(Maxwell)[p]*SigmaEII[p, i->j]         secondary-secondary EII
	//  eS[i][t] = Sum{j} * Sum{j} e[ij]*Int{_0^infty} dp*p^3*f(Maxwell)[p]*SigmaEII[p, i->j]  secondary-secondary EII energy loss
	// primary electron separate treatment:
	//  Gesc[t] = v[t] / R = sqrt(2 * Ep[t] / Np[t]) / R                                       primary electrons excaping beam zone
	//  Pht[i][t] = Sum{j} Sigma[i->j] * Intensity[t]
	//  ePht[i][t] = Sum{j} Sigma[i->j] * Intensity[t] * (omega - e[ij])
	//  Cll[t](N,Np,E,Ep)                                                                      primary-secondary energy exchange.
	//  W[i][t] = v[t] * Sum{j} Int{_0^(0.5*v^2 - EII.ionB[ij]) dw*w*dSigmaEII[v[t], i->j]/dw}
	//  Sp[i][t] = v[t] * Sum{j} SigmaEII[v[t], i->j]
	//  eSp[i][t] = v[t] * Sum{j} e[ij]*SigmaEII[v[t], i->j]
	adams_n = 5;
	t.resize(dt.size());
	if (f.size() != dt.size()) f = vector<double>(dt.size(), 0);

	int tot_p_size = 0;
	int p_size = adams_n + 1;
	for (auto& elem: Store) tot_p_size += elem.num_conf;

	dpdt.resize(tot_p_size);
	p.resize(tot_p_size);

	for (int i = 0; i < tot_p_size; i++) {
		dpdt[i].resize(p_size, 0.);
		p[i].resize(p_size, 0.);
	}

	int init_p = 0;
	for (auto& elem: Store) {
		p[init_p][0] = 1.;
		init_p += elem.num_conf;
	}

	vector<vector<double*>> map_p(Store.size());
	vector<vector<double*>> map_dpdt(Store.size());
	init_p = 0;
	for (int a = 0; a < Store.size(); a++) {
		for (int i = 0; i < Store[a].num_conf; i++) {
			map_p[a].push_back(p[i + init_p].data());
			map_dpdt[a].push_back(dpdt[i + init_p].data());
		}
		init_p += Store[a].num_conf;
	}

	// Initialize p, dpdt for the subsequent integration.
	init_p = 0;
	vector<vector<double>> A(Store.size());
	vector<vector<double>> X(Store.size());

	vector<vector<double>> Pht(Store.size());
	vector<vector<double>> ePht(Store.size());
	vector<vector<double>> Aug(Store.size());
	vector<vector<double>> eAug(Store.size());

	vector<vector<double>> S(Store.size());// Summed secondary EII rate.
	vector<vector<double>> eS(Store.size());// Summed secondary EII rate of electron energy loss (ion. potential).
	vector<vector<double>> Sp(Store.size());// Summed primary EII rate.
	vector<vector<double>> eSp(Store.size());// Summed EII rate of primary electron energy loss due to ionization potential.
	vector<vector<double>> W(Store.size());// Summed EII rate of primary electron energy loss to secondary electron creation.

	double tollerance = 0.000001, error = 1, old_p = 0, tmp = 0;
	double Temperature = 0, e_t = 0;

	for(int a = 0; a < Store.size(); a++) {
		A[a].resize(Store[a].num_conf, 0.);
		X[a].resize(Store[a].num_conf, 0.);
		Pht[a].resize(Store[a].num_conf, 0.);
		ePht[a].resize(Store[a].num_conf, 0.);
		Aug[a].resize(Store[a].num_conf, 0.);
		eAug[a].resize(Store[a].num_conf, 0.);
		S[a].resize(Store[a].num_conf, 0.);
		eS[a].resize(Store[a].num_conf, 0.);
		Sp[a].resize(Store[a].num_conf, 0.);
		eSp[a].resize(Store[a].num_conf, 0.);
		W[a].resize(Store[a].num_conf, 0.);

		for (auto& rate: Store[a].Auger) {
			Aug[a][rate.from] += rate.val;
			eAug[a][rate.from] += rate.val * rate.energy;
		}
		for (auto& rate: Store[a].Photo) {
			if (rate.from == 0) {
				e_t += rate.val * rate.energy;
				tmp += rate.val;
			}
		}
	}

	e_t /= tmp;
	double v_t = sqrt(2*e_t);
	double Gesc = 1.5*v_t/Store[0].R;
  double Pi4 = 4*Constant::Pi;
	tmp = 0;
	for (auto& at_Store: Store) {
		for (auto& eii: at_Store.EIIparams) {
			if (eii.init != 0) continue;
			for (int i = 0; i < eii.fin.size(); i++) {
				Temperature += Elecs.sigmaBEBw1(e_t, eii.ionB[i], eii.kin[i], eii.occ[i]) * at_Store.nAtoms;
				tmp += Elecs.sigmaBEB(e_t, eii.ionB[i], eii.kin[i], eii.occ[i]) * at_Store.nAtoms;
			}
			break;
		}
	}
	Temperature *= 2./3./tmp;
	double Factor = v_t;
	tmp = 0;
	Elecs.SetMaxwellPF(Temperature);
	double tmp_dNdt = 0, tmp_dEdt = 0, tmp_dNpdt = 0, tmp_dEpdt = 0;

	for (int m = 0; m < adams_n; m++) {
		if (m > 0) {
			for (auto& v : p) {// Guess.
				v[m] = v[m - 1];
			}

			Elecs.N[m] = Elecs.N[m-1];
			Elecs.E[m] = Elecs.E[m-1];
			Elecs.Np[m] = Elecs.Np[m-1];
			Elecs.Ep[m] = Elecs.Ep[m-1];
		}

		while (error > tollerance) {
			error = 0;

			Elecs.dNdt[m] = 0;
			Elecs.dEdt[m] = 0;
			Elecs.dNpdt[m] = 0;
			Elecs.dEpdt[m] = 0;

			for (int a = 0; a < Store.size(); a++) {			
		
				for (auto& rate: Store[a].Photo) {
					tmp = rate.val*f[m];
					A[a][rate.from] -= tmp;
					Pht[a][rate.from] += tmp;
					ePht[a][rate.from] += tmp * rate.energy;
					X[a][rate.to] += tmp * *(map_p[a][rate.from] + m);
				}
				for (auto& rate: Store[a].Fluor) {
					A[a][rate.from] -= rate.val;
					X[a][rate.to] += rate.val * *(map_p[a][rate.from] + m);
				}
				for (auto& rate: Store[a].Auger) {
					A[a][rate.from] -= rate.val;
					X[a][rate.to] += rate.val * *(map_p[a][rate.from] + m);
				}

				if (m > 0) {
					for (auto& EII: Store[a].EIIparams) {
						for (int i = 0; i < EII.fin.size(); i++) {
							tmp = Pi4*Elecs.MaxwellEII(EII.ionB[i], EII.kin[i], EII.occ[i]);
							S[a][EII.init] += tmp;
							eS[a][EII.init] += tmp*EII.ionB[i];
							old_p = Factor * Elecs.sigmaBEB(e_t, EII.ionB[i], EII.kin[i], EII.occ[i]);
							Sp[a][EII.init] += old_p;
							eSp[a][EII.init] += old_p * EII.ionB[i];
							X[a][EII.fin[i]] += *(map_p[a][EII.init] + m) * ( old_p * Elecs.Np[m] + tmp * Elecs.N[m] );
							W[a][EII.init] += Factor * Elecs.sigmaBEBw1(e_t, EII.ionB[i], EII.kin[i], EII.occ[i]);
						}				
					}
				}

				tmp_dEdt = 0;
				tmp_dNdt = 0;
				tmp_dEpdt = 0;
				tmp_dNpdt = 0;
				
				// Correct atomic occupancies.
				for (int i = 0; i < Store[a].num_conf; i++)
				{
					A[a][i] -= S[a][i] * Elecs.N[m] + Sp[a][i] * Elecs.Np[m];
					if (A[a][i] == 0 && X[a][i] == 0) continue;
					*(map_dpdt[a][i] + m) = A[a][i]* *(map_p[a][i] + m) + X[a][i];

					A[a][i] = 0.;
					X[a][i] = 0.;

					if (m > 0)
					{
						old_p = *(map_p[a][i] + m);
						*(map_p[a][i] + m) = *(map_p[a][i] + m - 1) + dt[m - 1]*0.5*(*(map_dpdt[a][i] + m) + *(map_dpdt[a][i] + m - 1) );
						if (old_p == 0) tmp = fabs(*(map_p[a][i] + m) - old_p);
						else tmp = fabs(*(map_p[a][i] + m) - old_p) / *(map_p[a][i] + m);
						if (error < tmp) error = tmp;
					}

					tmp_dNdt += *(map_p[a][i] + m) * ( Aug[a][i] + Elecs.N[m]*S[a][i] + Elecs.Np[m]*Sp[a][i] );
					tmp_dEdt += *(map_p[a][i] + m) * ( eAug[a][i] - Elecs.N[m]*eS[a][i] + Elecs.Np[m]*W[a][i] );
					tmp_dNpdt += *(map_p[a][i] + m) * ( Pht[a][i] - Gesc*Elecs.Np[m]);
					tmp_dEpdt += *(map_p[a][i] + m) * ( ePht[a][i] - Gesc*Elecs.Ep[m] - Elecs.Np[m]*(eSp[a][i] + W[a][i]));

					Pht[a][i] = 0.;
					ePht[a][i] = 0.;

					S[a][i] = 0.;
					eS[a][i] = 0.;
					Sp[a][i] = 0.;
					eSp[a][i] = 0.;
					W[a][i] = 0.;		
				}

				// Correct plasma.
				Elecs.dNdt[m] += Store[a].nAtoms * tmp_dNdt;
				Elecs.dEdt[m] += Store[a].nAtoms * tmp_dEdt;
				Elecs.dNpdt[m] += Store[a].nAtoms * tmp_dNpdt;
				Elecs.dEpdt[m] += Store[a].nAtoms * tmp_dEpdt;
			}

			if (m > 0) {
				Elecs.N[m] = Elecs.N[m - 1] + dt[m - 1]*0.5*(Elecs.dNdt[m] + Elecs.dNdt[m-1]);
				Elecs.E[m] = Elecs.E[m - 1] + dt[m - 1]*0.5*(Elecs.dEdt[m] + Elecs.dEdt[m-1]);
				Elecs.Np[m] = Elecs.Np[m - 1] + dt[m - 1]*0.5*(Elecs.dNpdt[m] + Elecs.dNpdt[m-1]);
				Elecs.Ep[m] = Elecs.Ep[m - 1] + dt[m - 1]*0.5*(Elecs.dEpdt[m] + Elecs.dEpdt[m-1]);

				Factor = Temperature;
				if (Elecs.N[m] > 0) {
					Temperature = 2*Elecs.E[m]/3/Elecs.N[m];
					Elecs.SetMaxwellPF(Temperature);
					if (error < fabs(Temperature/Factor - 1) ) error = fabs(Temperature/Factor - 1);
				}

				e_t = Elecs.Ep[m]/Elecs.Np[m];
				v_t = sqrt(2*e_t);
				Factor = v_t * 0.25/Constant::Pi;
			}
		}
		error = 1;
	}
}


IntegrateRateEquation::IntegrateRateEquation(vector<double> &dT, vector<double> &T, AtomRateData& Store, Plasma & Elecs, vector<double> InitCond, const vector<double>& Intensity) : dt(dT), t(T), f(Intensity), store(Store)
{
	// f(F) is intensity defined at times T[m].
	// InitCond defines number of states and initial values for p.

	// Recast plasma equations (instant thermalization of photoelectrons)
	//
	//  dNp[t]/dt = n*Sum{i} p[i][t]*( Pht[i] - Gesc[t]*Np[t] )
	//  dEp[t]/dt = n*Sum{i} p[i][t]*( ePht[i] - Np[t]*(W[i] + eSp[i]) - Gesc[t]*Ep[t] - Cll[t](N,Np,E,Ep) )	
	//  dN[t]/dt = n*Sum{i} p[i][t]*( Aug[i][t] + N[t]*S[i][t] + Np[t]*Sp[i][t])
	//  dE[t]/dt = n*Sum{i} p[i][t]*( eAug[i][t] - N[t]*eS[i] + Np[t]*W[i])	
	//
	//	e[ij] - ionization potential of electron with ion going from "i" to "j"
	//  S[i][t] = Sum{j} * Sum{j} Int{_0^infty} dp*p^3*f(Maxwell)[p]*SigmaEII[p, i->j]         secondary-secondary EII
	//  eS[i][t] = Sum{j} * Sum{j} e[ij]*Int{_0^infty} dp*p^3*f(Maxwell)[p]*SigmaEII[p, i->j]  secondary-secondary EII energy loss
	// primary electron separate treatment:
	//  Gesc[t] = v[t] / R = sqrt(2 * Ep[t] / Np[t]) / R                                       primary electrons excaping beam zone
	//  Pht[i][t] = Sum{j} Sigma[i->j] * Intensity[t]
	//  ePht[i][t] = Sum{j} Sigma[i->j] * Intensity[t] * (omega - e[ij])
	//  Cll[t](N,Np,E,Ep)                                                                      primary-secondary energy exchange.
	//  W[i][t] = v[t] * Sum{j} Int{_0^(0.5*v^2 - EII.ionB[ij]) dw*w*dSigmaEII[v[t], i->j]/dw}
	//  Sp[i][t] = v[t] * Sum{j} SigmaEII[v[t], i->j]
	//  eSp[i][t] = v[t] * Sum{j} e[ij]*SigmaEII[v[t], i->j]
	adams_n = 5;
	t.resize(dt.size());
	if (f.size() != dt.size()) f = vector<double>(dt.size(), 0);
	dpdt.resize(InitCond.size());
	p.resize(InitCond.size());

	for (int i = 0; i < InitCond.size(); i++) {
		dpdt[i].resize(adams_n + 1);
		p[i].resize(adams_n + 1);
		p[i][0] = InitCond[i];
	}

	// Initialize p, dpdt for the subsequent integration.
	vector<double> A(p.size(), 0);
	vector<double> X(p.size(), 0);

	vector<double> Pht(p.size(), 0);
	vector<double> ePht(p.size(), 0);
	vector<double> Aug(p.size(), 0);
	vector<double> eAug(p.size(), 0);

	vector<double> S(p.size(), 0);// Summed secondary EII rate.
	vector<double> eS(p.size(), 0);// Summed secondary EII rate of electron energy loss (ion. potential).
	vector<double> Sp(p.size(), 0);// Summed primary EII rate.
	vector<double> eSp(p.size(), 0);// Summed EII rate of primary electron energy loss due to ionization potential.
	vector<double> W(p.size(), 0);// Summed EII rate of primary electron energy loss to secondary electron creation.

	double tollerance = 0.000001, error = 1, old_p = 0, tmp = 0;
	double Temperature = 0;

	for (auto& rate: Store.Auger) {
			Aug[rate.from] += rate.val;
			eAug[rate.from] += rate.val * rate.energy;
	}
	double e_t = 0;
	for (auto& rate: Store.Photo) {
		if (rate.from == 0) {
			e_t += rate.val * rate.energy;
			tmp += rate.val;
		}
	}
	e_t /= tmp;
	double v_t = sqrt(2*e_t);
	double Gesc = v_t/Store.R;
	tmp = 0;
	for (auto& eii: Store.EIIparams) {
		if (eii.init != 0) continue;
		for (int i = 0; i < eii.fin.size(); i++) {
			Temperature += Elecs.sigmaBEBw1(e_t, eii.ionB[i], eii.kin[i], eii.occ[i]);
			tmp += Elecs.sigmaBEB(e_t, eii.ionB[i], eii.kin[i], eii.occ[i]);
		}
		break;
	}
	Temperature *= 2./3./tmp;
	double Factor = v_t * 0.25/Constant::Pi;
	tmp = 0;
	Elecs.SetMaxwellPF(Temperature);

	for (int m = 0; m < adams_n; m++) {
		if (m > 0) {
			for (auto& v : p) {// Guess.
				v[m] = v[m - 1];
			}

			Elecs.N[m] = Elecs.N[m-1];
			Elecs.E[m] = Elecs.E[m-1];
			Elecs.Np[m] = Elecs.Np[m-1];
			Elecs.Ep[m] = Elecs.Ep[m-1];
		}

		while (error > tollerance) {
			error = 0;

			Elecs.dNdt[m] = 0;
			Elecs.dEdt[m] = 0;
			Elecs.dNpdt[m] = 0;
			Elecs.dEpdt[m] = 0;

			for (auto& rate: Store.Photo) {
				tmp = rate.val*f[m];
				A[rate.from] -= tmp;
				Pht[rate.from] += tmp;
				ePht[rate.from] += tmp * rate.energy;
				X[rate.to] += tmp * p[rate.from][m];
			}
			for (auto& rate: Store.Fluor) {
				A[rate.from] -= rate.val;
				X[rate.to] += rate.val * p[rate.from][m];
			}
			for (auto& rate: Store.Auger) {
				A[rate.from] -= rate.val;
				X[rate.to] += rate.val * p[rate.from][m];
			}

			if (m > 0) {
				for (auto& EII: Store.EIIparams) {
					for (int i = 0; i < EII.fin.size(); i++) {
						tmp = Elecs.MaxwellEII(EII.ionB[i], EII.kin[i], EII.occ[i]);
						S[EII.init] += tmp;
						eS[EII.init] += tmp*EII.ionB[i];
						old_p = Factor * Elecs.sigmaBEB(e_t, EII.ionB[i], EII.kin[i], EII.occ[i]);
						Sp[EII.init] += old_p;
						eSp[EII.init] += old_p * EII.ionB[i];
						X[EII.fin[i]] += p[EII.init][m] * ( old_p * Elecs.Np[m] + tmp * Elecs.N[m] );
						W[EII.init] += Factor * Elecs.sigmaBEBw1(e_t, EII.ionB[i], EII.kin[i], EII.occ[i]);
					}				
				}
			}

			for (int i = 0; i < p.size(); i++)// Correct atomic occupancies.
			{
				A[i] -= S[i] * Elecs.N[m] + Sp[i] * Elecs.Np[m];
				if (A[i] == 0 && X[i] == 0) continue;
				dpdt[i][m] = A[i]*p[i][m] + X[i];

				A[i] = 0.;
				X[i] = 0.;

				if (m > 0)
				{
					old_p = p[i][m];
					p[i][m] = p[i][m - 1] + dt[m - 1]*0.5*(dpdt[i][m] + dpdt[i][m-1]);
					if (old_p == 0) tmp = fabs(p[i][m] - old_p);
					else tmp = fabs(p[i][m] - old_p) / p[i][m];
					if (error < tmp) error = tmp;
				}

				Elecs.dNdt[m] += p[i][m] * ( Aug[i] + Elecs.N[m]*S[i] + Elecs.Np[m]*Sp[i] );
				Elecs.dEdt[m] += p[i][m] * ( eAug[i] - Elecs.N[m]*eS[i] + Elecs.Np[m]*W[i] );
				Elecs.dNpdt[m] += p[i][m] * ( Pht[i] - Gesc*Elecs.Np[m]);
				Elecs.dEpdt[m] += p[i][m] * ( ePht[i] - Gesc*Elecs.Ep[m] - Elecs.Np[m]*(eSp[i] + W[i]));

				Pht[i] = 0.;
				ePht[i] = 0.;

				S[i] = 0.;
				eS[i] = 0.;
				Sp[i] = 0.;
				eSp[i] = 0.;
				W[i] = 0.;		
			}// Correct plasma.
			Elecs.dNdt[m] *= store.nAtoms;
			Elecs.dEdt[m] *= store.nAtoms;
			Elecs.dNpdt[m] *= store.nAtoms;
			Elecs.dEpdt[m] *= store.nAtoms;
			if (m > 0) {
				Elecs.N[m] = Elecs.N[m - 1] + dt[m - 1]*0.5*(Elecs.dNdt[m] + Elecs.dNdt[m-1]);
				Elecs.E[m] = Elecs.E[m - 1] + dt[m - 1]*0.5*(Elecs.dEdt[m] + Elecs.dEdt[m-1]);
				Elecs.Np[m] = Elecs.Np[m - 1] + dt[m - 1]*0.5*(Elecs.dNpdt[m] + Elecs.dNpdt[m-1]);
				Elecs.Ep[m] = Elecs.Ep[m - 1] + dt[m - 1]*0.5*(Elecs.dEpdt[m] + Elecs.dEpdt[m-1]);
				e_t = Elecs.Ep[m]/Elecs.Np[m];
				v_t = sqrt(2*e_t);
				Factor = v_t * 0.25/Constant::Pi;

				if (Elecs.N[m] > 0) {
					Temperature = 2*Elecs.E[m]/3/Elecs.N[m];
					Elecs.SetMaxwellPF(Temperature);
				}
			}
		}
		error = 1;
	}
}


int IntegrateRateEquation::Solve(double P_min, double P_max, int storage_time_pts)
{
	// Set up storage container for time (t_storage) and occupancies (p_storage)
	int storage_count = t.size() / storage_time_pts;
	time_storage.resize(t.size()/storage_count);
	p_storage.resize(p.size());
	for (auto& v : p_storage) {
		v.resize(t.size() / storage_count);
	}

	for (int m = 0; m < adams_n; m++) {
		if ((m % storage_count) == 0) {
			int k = m / storage_count;
			time_storage[k] = t[m];
			for (int i = 0; i < p.size(); i++) {
				p_storage[i][k] = p[i][m];
			}
		}
	}

	//Begin numerical integration with classic Adams-Moulton predictor-corrector scheme
	vector<double> A(p.size()), X(p.size());
	double tmp = 0;
	bool unstable = false;
	int p_size = p[0].size();
	for (int m = adams_n; m < t.size(); m++) {
		cout << "\r" << m << "/" << t.size();
		/*for (int i = 0; i < p.size(); i++) {
			p[i].push_back(0);
			dpdt[i].push_back(0);
		}*/
		for (vector<double>::iterator it = A.begin(); it != A.end(); it++) *it = 0;
		for (vector<double>::iterator it = X.begin(); it != X.end(); it++) *it = 0;

		for (int i = 0; i < p.size(); i++) { // Predict p.
			tmp = p[i][adams_n-1];
			for (int j = 0; j < adams_n; j++) {
				tmp += Bashforth_5[j] * dpdt[i][adams_n - j - 1] * dt[m - j - 1];
			}
			p[i].back() = tmp;
		}

		for (auto& rate: store.Photo) {//predict dpdt
			tmp = rate.val*f[m];
			A[rate.from] -= tmp;
			X[rate.to] += tmp * p[rate.from].back();
		}
		for (auto& rate: store.Fluor) {
			A[rate.from] -= rate.val;
			X[rate.to] += rate.val * p[rate.from].back();
		}
		for (auto& rate: store.Auger) {
			A[rate.from] -= rate.val;
			X[rate.to] += rate.val * p[rate.from].back();
		}
		/*
		for (auto& rate: store.Rates) {
			if (CompareChar(depends_t, rate.type)) {
				A[rate.from] -= rate.val*f[m];
				X[rate.to] += rate.val * f[m] * p[rate.from].back();
			} else {
				A[rate.from] -= rate.val;
				X[rate.to] += rate.val * p[rate.from].back();
			}
		}
		*/
		for (int i = 0; i < p.size(); i++) {// Free up last element in p[i] and erase the first one.
			dpdt[i].back() = A[i] * p[i].back() + X[i];
			for (int j = 1; j < p_size; j++) {
				p[i][j - 1] = p[i][j];
				dpdt[i][j - 1] = dpdt[i][j];
			}
			p[i].back() = 0;
			dpdt[i].back() = 0;
		}
		for (vector<double>::iterator it = X.begin(); it != X.end(); it++) *it = 0;

		for (int i = 0; i < p.size(); i++) {// Correct p.
			tmp = p[i][adams_n-2];
			for (int j = 0; j < adams_n; j++) {
				tmp += Moulton_5[j] * dpdt[i][adams_n - j - 1] * dt[m - j];
			}
			p[i][p_size - 2] = tmp;
			if (tmp < P_min || tmp > P_max) unstable = true;
		}

		for (auto& rate: store.Photo) {//correct dpdt
			X[rate.to] += rate.val * f[m] * p[rate.from][p_size - 2];
		}
		for (auto& rate: store.Fluor) {
			X[rate.to] += rate.val * p[rate.from][p_size - 2];
		}
		for (auto& rate: store.Auger) {
			X[rate.to] += rate.val * p[rate.from][p_size - 2];
		}

		if (unstable) {
			cout << "\r";
			return m;
		}

		for (int i = 0; i < p.size(); i++) {
			dpdt[i][p_size - 2] = A[i] * p[i][p_size - 2] + X[i];
		}

		if ((m % storage_count) == 0) {
			int k = m / storage_count;
			time_storage[k] = t[m];
			for (int i = 0; i < p.size(); i++)	{
				p_storage[i][k] = p[i][p_size - 2];
			}
		}
	}

	cout << endl;
	return 0;
}

int IntegrateRateEquation::Solve(Plasma & Elecs, double P_min, double P_max, int storage_time_pts)
{
	// f(F) is intensity defined at times T[m].
	// InitCond defines number of states and initial values for p.

	// Recast plasma equations (instant thermalization of photoelectrons)
	//
	//  dNp[t]/dt = n*Sum{i} p[i][t]*( Pht[i] - Gesc[t]*Np[t] )
	//  dEp[t]/dt = n*Sum{i} p[i][t]*( ePht[i] - Np[t]*(W[i] + eSp[i]) - Gesc[t]*Ep[t] - Cll[t](N,Np,E,Ep) )	
	//  dN[t]/dt = n*Sum{i} p[i][t]*( Aug[i][t] + N[t]*S[i][t] + Np[t]*Sp[i][t] - N[t]*N[t]*Tbr[i][t] )
	//  dE[t]/dt = n*Sum{i} p[i][t]*( eAug[i][t] - N[t]*eS[i] + Np[t]*W[i] + N[t]*N[t]*eTbr[i][t])	
	//
	// e[ij] - ionization potential of electron with ion going from "i" to "j"
	//  S[i][t] = Sum{j} Int{_0^infty} dp*p^3*f(Maxwell)[p]*SigmaEII[p, i->j] = Sum{j} S[ij][t] ( "t" dependence through temperature)
	//  eS[i][t] = Sum{j} e[ij]*Int{_0^infty} dp*p^3*f(Maxwell)[p]*SigmaEII[p, i->j]
	//  Tbr[i][t] = Sum{j} (2\pi)^3*exp(e[ij]/T)*S[ij][t] - three body recombination
	//  eTbr[i][t] = Sum{j} (2\pi)^3*exp(e[ij]/T)*eS[ij][t]
	// primary electron separate treatment:
	//  Gesc[t] = v[t] / R = sqrt(2 * Ep[t] / Np[t]) / R                                       primary electrons excaping beam zone
	//  Pht[i][t] = Sum{j} Sigma[i->j] * Intensity[t]
	//  ePht[i][t] = Sum{j} Sigma[i->j] * Intensity[t] * (omega - e[ij])
	//  Cll[t](N,Np,E,Ep)                                                                      primary-secondary energy exchange.
	//  W[i][t] = v[t] * Sum{j} Int{_0^(0.5*v^2 - EII.ionB[ij]) dw*w*dSigmaEII[v[t], i->j]/dw}
	//  Sp[i][t] = v[t] * Sum{j} SigmaEII[v[t], i->j]
	//  eSp[i][t] = v[t] * Sum{j} e[ij]*SigmaEII[v[t], i->j]

	// Set up storage container for time (t_storage) and occupancies (p_storage)
	int storage_count = t.size() / storage_time_pts;
	time_storage.resize(t.size()/storage_count);
	p_storage.resize(p.size());
	for (auto& v : p_storage) {
		v.resize(t.size() / storage_count);
	}

	for (int m = 0; m < adams_n; m++) {
		if ((m % storage_count) == 0) {
			int k = m / storage_count;
			time_storage[k] = t[m];
			for (int i = 0; i < p.size(); i++) {
				p_storage[i][k] = p[i][m];
			}
		}
	}

	//Begin numerical integration with classic Adams-Moulton predictor-corrector scheme
	vector<double> A(p.size(), 0.), X(p.size(), 0.);

	vector<double> Pht(p.size(), 0.);
	vector<double> ePht(p.size(), 0.);
	vector<double> Aug(p.size(), 0.);
	vector<double> eAug(p.size(), 0.);

	vector<double> S(p.size(), 0.);// Summed EII rate.
	vector<double> eS(p.size(), 0.);// Summed EII rate of electron energy loss.
	vector<double> Tbr(p.size(), 0.);// Summed TBR rate.
	vector<double> eTbr(p.size(), 0.);// Summed TBR rate of electron energy loss.
	vector<double> Sp(p.size(), 0);// Summed photo EII rate.
	vector<double> eSp(p.size(), 0);// Summed EII rate of photo-electron energy loss due to ionization potential.
	vector<double> W(p.size(), 0);// Summed EII rate of photo electron energy loss to secondary electron creation.

	double e_t = Elecs.Ep[adams_n-1]/Elecs.Np[adams_n-1];
	double v_t = sqrt(2*e_t);
	double Factor = v_t * 0.25/Constant::Pi;
	double Gesc = v_t/store.R;
	double Temperature = 2*Elecs.E[adams_n-1]/3/Elecs.N[adams_n-1];
	Elecs.SetMaxwellPF(Temperature);

	for (auto& rate: store.Auger) {
		Aug[rate.from] += rate.val;
		eAug[rate.from] += rate.val * rate.energy;
	}

	double tmp = 0, mxw_tmp = 0, Pi3x8 = 8*pow(Constant::Pi, 3);
	int p_size = p[0].size();

	for (int m = adams_n; m < t.size(); m++) {
		cout << "\r" << m << "/" << t.size();
		
		for (int i = 0; i < A.size(); i++) {
			A[i] = 0.;
			X[i] = 0.;
			Pht[i] = 0.;
			ePht[i] = 0.;
			S[i] = 0.;
			eS[i] = 0.;
			Tbr[i] = 0.;
			eTbr[i] = 0.;
			Sp[i] = 0.;
			eSp[i] = 0;
			W[i] = 0.;
		}

		for (int i = 0; i < p.size(); i++) { // Predict p.
			tmp = p[i][adams_n-1];
			for (int j = 0; j < adams_n; j++) {
				tmp += Bashforth_5[j] * dpdt[i][adams_n - j - 1] * dt[m - j - 1];
			}
			p[i].back() = tmp;
		}
		Elecs.E[m] = Elecs.E[m-1];// Predict plasma E and N.
		Elecs.N[m] = Elecs.N[m-1];
		Elecs.Ep[m] = Elecs.Ep[m-1];// Predict plasma E and N.
		Elecs.Np[m] = Elecs.Np[m-1];
		for (int j = 0; j < adams_n; j++) {
			Elecs.E[m] += Bashforth_5[j] * Elecs.dEdt[m - j - 1] * dt[m - j - 1];
			Elecs.N[m] += Bashforth_5[j] * Elecs.dNdt[m - j - 1] * dt[m - j - 1];
			Elecs.Ep[m] += Bashforth_5[j] * Elecs.dEpdt[m - j - 1] * dt[m - j - 1];
			Elecs.Np[m] += Bashforth_5[j] * Elecs.dNpdt[m - j - 1] * dt[m - j - 1];
		}
		e_t = Elecs.Ep[m] / Elecs.Np[m];
		v_t = sqrt(2*e_t);
		Gesc = v_t / store.R;
		Factor = v_t * 0.25 / Constant::Pi;	
		Temperature = 2*Elecs.E[m]/3/Elecs.N[m];
		Elecs.SetMaxwellPF(Temperature);

		for (auto& rate: store.Photo) {//predict dpdt
			tmp = rate.val * f[m];
			A[rate.from] -= tmp;
			Pht[rate.from] += tmp;
			ePht[rate.from] += tmp * rate.energy;
			X[rate.to] += tmp * p[rate.from].back();
		}
		for (auto& rate: store.Fluor) {
			A[rate.from] -= rate.val;
			X[rate.to] += rate.val * p[rate.from].back();
		}
		for (auto& rate: store.Auger) {
			A[rate.from] -= rate.val;
			X[rate.to] += rate.val * p[rate.from].back();
		}
		for (auto& EII: store.EIIparams) {
			for (int i = 0; i < EII.fin.size(); i++) {
				tmp = Elecs.MaxwellEII(EII.ionB[i], EII.kin[i], EII.occ[i]);
				mxw_tmp = Pi3x8 * Elecs.MaxwellPF(EII.ionB[i]) * tmp;
				Tbr[EII.fin[i]] += mxw_tmp;
				eTbr[EII.fin[i]] += mxw_tmp * EII.ionB[i];
				X[EII.init] += mxw_tmp * p[EII.fin[i]].back() * Elecs.N[m] * Elecs.N[m];
				S[EII.init] += tmp;
				eS[EII.init] += tmp * EII.ionB[i];

				mxw_tmp = Factor * Elecs.sigmaBEB(e_t, EII.ionB[i], EII.kin[i], EII.occ[i]);
				Sp[EII.init] += mxw_tmp;
				eSp[EII.init] += mxw_tmp * EII.ionB[i];
				X[EII.fin[i]] += p[EII.init].back() * (mxw_tmp * Elecs.Np[m] + tmp * Elecs.N[m]);
				W[EII.init] += Factor * Elecs.sigmaBEBw1(e_t, EII.ionB[i], EII.kin[i], EII.occ[i]);
			}		
		}
		
		for (int i = 0; i < p.size(); i++) {// Free up last element in p[i] and erase the first one.
			A[i] -= Elecs.N[m] * (S[i] + Elecs.N[m] * Tbr[i]) + Elecs.Np[m] * Sp[i];// 
			tmp = A[i] * p[i].back() + X[i];
			dpdt[i].back() = tmp;
			Elecs.dNpdt[m] += p[i].back() * ( Pht[i] - Gesc * Elecs.Np[m] );
			Elecs.dEpdt[m] += p[i].back() * ( ePht[i] - Gesc * Elecs.Ep[m] - Elecs.Np[m] * (eSp[i] + W[i]) );
			Elecs.dNdt[m] += p[i].back() * ( Aug[i] + Elecs.N[m] * (S[i] - Elecs.N[m] * Tbr[i]) + Elecs.Np[m] * Sp[i] );
			Elecs.dEdt[m] += p[i].back() * ( eAug[i] - Elecs.N[m] * (eS[i] - Elecs.N[m] * eTbr[i]) + Elecs.Np[m] * W[i] );

			for (int j = 1; j < p_size; j++) {
				p[i][j - 1] = p[i][j];
				dpdt[i][j - 1] = dpdt[i][j];
			}
			p[i].back() = 0;
			dpdt[i].back() = 0;
		}
		Elecs.dNdt[m] *= store.nAtoms;
		Elecs.dEdt[m] *= store.nAtoms;
		Elecs.dNpdt[m] *= store.nAtoms;
		Elecs.dEpdt[m] *= store.nAtoms;

		for (int i = 0; i < p.size(); i++) {// Correct p.
			tmp = p[i][adams_n - 2];
			for (int j = 0; j < adams_n; j++) {
				tmp += Moulton_5[j] * dpdt[i][adams_n - j - 1] * dt[m - j];
			}
			p[i][adams_n - 1] = tmp;
			if (tmp < P_min || tmp > P_max) {
				printf("\n i = %d \n", i);
				return m;
			}
		}

		for (int i = 0; i < A.size(); i++) {
			A[i] = 0.;
			X[i] = 0.;
			Pht[i] = 0.;
			ePht[i] = 0.;
			S[i] = 0.;
			eS[i] = 0.;
			Tbr[i] = 0.;
			eTbr[i] = 0.;
			Sp[i] = 0.;
			eSp[i] = 0;
			W[i] = 0.;
		}

		Elecs.E[m] = Elecs.E[m-1];// Correct plasma.
		Elecs.N[m] = Elecs.N[m-1];
		Elecs.Ep[m] = Elecs.Ep[m-1];
		Elecs.Np[m] = Elecs.Np[m-1];
		for (int j = 0; j < adams_n; j++) {
			Elecs.E[m] += Moulton_5[j] * Elecs.dEdt[m - j] * dt[m - j];
			Elecs.N[m] += Moulton_5[j] * Elecs.dNdt[m - j] * dt[m - j];
			Elecs.Ep[m] += Moulton_5[j] * Elecs.dEpdt[m - j] * dt[m - j];
			Elecs.Np[m] += Moulton_5[j] * Elecs.dNpdt[m - j] * dt[m - j];
		}
		e_t = Elecs.Ep[m] / Elecs.Np[m];
		v_t = sqrt(2*e_t);
		Gesc = v_t / store.R;
		Factor = v_t * 0.25 / Constant::Pi;
		Temperature = 2*Elecs.E[m]/3/Elecs.N[m];
		Elecs.SetMaxwellPF(Temperature);

		for (auto& EII: store.EIIparams) {
			for (int i = 0; i < EII.fin.size(); i++) {
				tmp = Elecs.MaxwellEII(EII.ionB[i], EII.kin[i], EII.occ[i]);
				mxw_tmp = Pi3x8 * Elecs.MaxwellPF(EII.ionB[i]) * tmp;
				Tbr[EII.fin[i]] += mxw_tmp;
				eTbr[EII.fin[i]] += mxw_tmp * EII.ionB[i];
				X[EII.init] += mxw_tmp * p[EII.fin[i]][p_size - 2] * Elecs.N[m] * Elecs.N[m];
				S[EII.init] += tmp;
				eS[EII.init] += tmp * EII.ionB[i];

				mxw_tmp = Factor * Elecs.sigmaBEB(e_t, EII.ionB[i], EII.kin[i], EII.occ[i]);
				Sp[EII.init] += mxw_tmp;
				eSp[EII.init] += mxw_tmp * EII.ionB[i];
				X[EII.fin[i]] += p[EII.init][p_size - 2] * ( mxw_tmp * Elecs.Np[m] + tmp * Elecs.N[m] );
				W[EII.init] += Factor * Elecs.sigmaBEBw1(e_t, EII.ionB[i], EII.kin[i], EII.occ[i]);
			}		
		}
		
		for (auto& rate: store.Photo) {//correct dpdt
			tmp = rate.val*f[m];
			A[rate.from] -= tmp;
			Pht[rate.from] += tmp;
			ePht[rate.from] += tmp * rate.energy;
			X[rate.to] += tmp * p[rate.from][p_size - 2];
		}
		for (auto& rate: store.Fluor) {
			A[rate.from] -= rate.val;
			X[rate.to] += rate.val * p[rate.from][p_size - 2];
		}
		for (auto& rate: store.Auger) {
			A[rate.from] -= rate.val;
			X[rate.to] += rate.val * p[rate.from][p_size - 2];
		}
		
		Elecs.dNdt[m] = 0.;// Correct plasma.
		Elecs.dEdt[m] = 0.;
		Elecs.dNpdt[m] = 0.;
		Elecs.dEpdt[m] = 0.;

		tmp = 0.;

		for (int i = 0; i < p.size(); i++) {
			A[i] -= Elecs.N[m] * (S[i] + Elecs.N[m] * Tbr[i]) + Elecs.Np[m] * Sp[i];
			dpdt[i][p_size - 2] = A[i] * p[i][p_size - 2] + X[i];
			Elecs.dNpdt[m] += p[i][p_size - 2] * ( Pht[i] - Gesc * Elecs.Np[m] );
			Elecs.dEpdt[m] += p[i][p_size - 2] * ( ePht[i] - Gesc * Elecs.Ep[m] - Elecs.Np[m] * (eSp[i] + W[i]) );
			Elecs.dNdt[m] += p[i][p_size - 2] * ( Aug[i] + Elecs.N[m] * (S[i] - Elecs.N[m] * Tbr[i]) + Elecs.Np[m] * Sp[i] );
			Elecs.dEdt[m] += p[i][p_size - 2] * ( eAug[i] - Elecs.N[m] * (eS[i] - Elecs.N[m] * eTbr[i]) + Elecs.Np[m] * W[i] );				
		}
		Elecs.dNdt[m] *= store.nAtoms;
		Elecs.dEdt[m] *= store.nAtoms;
		Elecs.dNpdt[m] *= store.nAtoms;
		Elecs.dEpdt[m] *= store.nAtoms;

		if ((m % storage_count) == 0) {
			int k = m / storage_count;
			time_storage[k] = t[m];
			for (int i = 0; i < p.size(); i++)	{
				p_storage[i][k] = p[i][p_size - 2];
			}
		}
	}

	cout << endl;
	return 0;
}

int IntegrateRateEquation::Solve(Plasma & Elecs, vector<AtomRateData> & Store, int storage_time_pts)
{
	// f(F) is intensity defined at times T[m].
	// InitCond defines number of states and initial values for p.

	// Recast plasma equations (instant thermalization of photoelectrons)
	//
	//  dNp[t]/dt = n*Sum{i} p[i][t]*( Pht[i] - Gesc[t]*Np[t] )
	//  dEp[t]/dt = n*Sum{i} p[i][t]*( ePht[i] - Np[t]*(W[i] + eSp[i]) - Gesc[t]*Ep[t] - Cll[t](N,Np,E,Ep) )	
	//  dN[t]/dt = n*Sum{i} p[i][t]*( Aug[i][t] + N[t]*S[i][t] + Np[t]*Sp[i][t] - N[t]*N[t]*Tbr[i][t] )
	//  dE[t]/dt = n*Sum{i} p[i][t]*( eAug[i][t] - N[t]*eS[i] + Np[t]*W[i] + N[t]*N[t]*eTbr[i][t])	
	//
	// e[ij] - ionization potential of electron with ion going from "i" to "j"
	//  S[i][t] = 4pi*Sum{j} Int{_0^infty} dp*p^3*f(Maxwell)[p]*SigmaEII[p, i->j] = Sum{j} S[ij][t] ( "t" dependence through temperature)
	//  eS[i][t] = 4pi*Sum{j} e[ij]*Int{_0^infty} dp*p^3*f(Maxwell)[p]*SigmaEII[p, i->j]
	//  Tbr[i][t] = Sum{j} (2\pi)^3*exp(e[ij]/T)*S[ij][t] - three body recombination
	//  eTbr[i][t] = Sum{j} (2\pi)^3*exp(e[ij]/T)*eS[ij][t]
	// primary electron separate treatment:
	//  Gesc[t] = v[t] / R = sqrt(2 * Ep[t] / Np[t]) / R                                       primary electrons excaping beam zone
	//  Pht[i][t] = Sum{j} Sigma[i->j] * Intensity[t]
	//  ePht[i][t] = Sum{j} Sigma[i->j] * Intensity[t] * (omega - e[ij])
	//  Cll[t](N,Np,E,Ep)                                                                      primary-secondary energy exchange.
	//  W[i][t] = v[t] * Sum{j} Int{_0^(0.5*v^2 - EII.ionB[ij]) dw*w*dSigmaEII[v[t], i->j]/dw}
	//  Sp[i][t] = v[t] * Sum{j} SigmaEII[v[t], i->j]
	//  eSp[i][t] = v[t] * Sum{j} e[ij]*SigmaEII[v[t], i->j]

	// Set up storage container for time (t_storage) and occupancies (p_storage)
	int storage_count = t.size() / storage_time_pts;
	time_storage.resize(t.size()/storage_count);
	p_storage.resize(p.size());
	for (auto& v : p_storage) {
		v.resize(t.size() / storage_count);
	}

	for (int m = 0; m < adams_n; m++) {
		if ((m % storage_count) == 0) {
			int k = m / storage_count;
			time_storage[k] = t[m];
			for (int i = 0; i < p.size(); i++) {
				p_storage[i][k] = p[i][m];
			}
		}
	}

	vector<vector<double*>> map_p(Store.size());
	vector<vector<double*>> map_dpdt(Store.size());
	int init_p = 0;
	for (int a = 0; a < Store.size(); a++) {
		for (int i = 0; i < Store[a].num_conf; i++) {
			map_p[a].push_back(p[i + init_p].data());
			map_dpdt[a].push_back(dpdt[i + init_p].data());
		}
		init_p += Store[a].num_conf;
	}

	// Initialize p, dpdt for the subsequent integration.
	init_p = 0;
	vector<vector<double>> A(Store.size());
	vector<vector<double>> X(Store.size());

	vector<vector<double>> Pht(Store.size());
	vector<vector<double>> ePht(Store.size());
	vector<vector<double>> Aug(Store.size());
	vector<vector<double>> eAug(Store.size());

	vector<vector<double>> S(Store.size());// Summed secondary EII rate.
	vector<vector<double>> eS(Store.size());// Summed secondary EII rate of electron energy loss (ion. potential).
	vector<vector<double>> Sp(Store.size());// Summed primary EII rate.
	vector<vector<double>> eSp(Store.size());// Summed EII rate of primary electron energy loss due to ionization potential.
	vector<vector<double>> W(Store.size());// Summed EII rate of primary electron energy loss to secondary electron creation.
	vector<vector<double>> Tbr(Store.size());//Three-body recombination.
	vector<vector<double>> eTbr(Store.size());// Energy transport associated with Tbr.

	double e_t = Elecs.Ep[adams_n-1]/Elecs.Np[adams_n-1];
	double v_t = sqrt(2*e_t);
	double Gesc = 1.*v_t/store.R;
	double Temperature = 2*Elecs.E[adams_n-1]/3/Elecs.N[adams_n-1];
	Elecs.SetMaxwellPF(Temperature);

	for(int a = 0; a < Store.size(); a++) {
		A[a].resize(Store[a].num_conf, 0.);
		X[a].resize(Store[a].num_conf, 0.);
		Pht[a].resize(Store[a].num_conf, 0.);
		ePht[a].resize(Store[a].num_conf, 0.);
		Aug[a].resize(Store[a].num_conf, 0.);
		eAug[a].resize(Store[a].num_conf, 0.);
		S[a].resize(Store[a].num_conf, 0.);
		eS[a].resize(Store[a].num_conf, 0.);
		Sp[a].resize(Store[a].num_conf, 0.);
		eSp[a].resize(Store[a].num_conf, 0.);
		W[a].resize(Store[a].num_conf, 0.);
		Tbr[a].resize(Store[a].num_conf, 0.);
		eTbr[a].resize(Store[a].num_conf, 0.);

		for (auto& rate: Store[a].Auger) {
			Aug[a][rate.from] += rate.val;
			eAug[a][rate.from] += rate.val * rate.energy;
		}
	}

	double tmp = 0, mxw_tmp = 0, Pi3x8 = 8*pow(Constant::Pi, 3), Pix4 = 4*Constant::Pi;
	int p_size = p[0].size();
	int tot_p_size = p.size();
	double tmp_dNdt = 0, tmp_dEdt = 0, tmp_dNpdt = 0, tmp_dEpdt = 0;

	for (int m = adams_n; m < t.size(); m++) {
		cout << "\r" << m << "/" << t.size();

		// Predict plasma E and N.
		Elecs.E[m] = Elecs.E[m-1];
		Elecs.N[m] = Elecs.N[m-1];
		Elecs.Ep[m] = Elecs.Ep[m-1];
		Elecs.Np[m] = Elecs.Np[m-1];
		for (int j = 0; j < adams_n; j++) {
			Elecs.E[m] += Bashforth_5[j] * Elecs.dEdt[m - j - 1] * dt[m - j - 1];
			Elecs.N[m] += Bashforth_5[j] * Elecs.dNdt[m - j - 1] * dt[m - j - 1];
			Elecs.Ep[m] += Bashforth_5[j] * Elecs.dEpdt[m - j - 1] * dt[m - j - 1];
			Elecs.Np[m] += Bashforth_5[j] * Elecs.dNpdt[m - j - 1] * dt[m - j - 1];
		}
		e_t = Elecs.Ep[m] / Elecs.Np[m];
		v_t = sqrt(2*e_t);
		Gesc = 1.*v_t / store.R;
		Temperature = 2*Elecs.E[m]/3/Elecs.N[m];
		Elecs.SetMaxwellPF(Temperature);

		for (int a = 0; a < Store.size(); a++) {
			tot_p_size = Store[a].num_conf;

			for (int i = 0; i < tot_p_size; i++) {
				A[a][i] = 0.;
				X[a][i] = 0.;
				Pht[a][i] = 0.;
				ePht[a][i] = 0.;
				S[a][i] = 0.;
				eS[a][i] = 0.;
				Tbr[a][i] = 0.;
				eTbr[a][i] = 0.;
				Sp[a][i] = 0.;
				eSp[a][i] = 0.;
				W[a][i] = 0.;

				// Predict p.
				tmp = *(map_p[a][i] + adams_n - 1);
				for (int j = 0; j < adams_n; j++) {
					tmp += Bashforth_5[j] * *(map_dpdt[a][i] + adams_n - j - 1) * dt[m - j - 1];
				}
				*(map_p[a][i] + adams_n) = tmp;
			}
			//predict dpdt
			for (auto& rate: Store[a].Photo) {
				tmp = rate.val * f[m];
				A[a][rate.from] -= tmp;
				Pht[a][rate.from] += tmp;
				ePht[a][rate.from] += tmp * rate.energy;
				X[a][rate.to] += tmp * *(map_p[a][rate.from] + adams_n);
			}
			for (auto& rate: Store[a].Fluor) {
				A[a][rate.from] -= rate.val;
				X[a][rate.to] += rate.val * *(map_p[a][rate.from] + adams_n);
			}
			for (auto& rate: Store[a].Auger) {
				A[a][rate.from] -= rate.val;
				X[a][rate.to] += rate.val * *(map_p[a][rate.from] + adams_n);
			}
			for (auto& EII: Store[a].EIIparams) {
				for (int i = 0; i < EII.fin.size(); i++) {
					tmp = Pix4*Elecs.MaxwellEII(EII.ionB[i], EII.kin[i], EII.occ[i]);
					mxw_tmp = Pi3x8 * Elecs.MaxwellPF(-1.*EII.ionB[i]) * tmp;
					Tbr[a][EII.fin[i]] += mxw_tmp;
					eTbr[a][EII.fin[i]] += mxw_tmp * EII.ionB[i];
					X[a][EII.init] += mxw_tmp * *(map_p[a][EII.fin[i]] + adams_n) * Elecs.N[m] * Elecs.N[m];
					S[a][EII.init] += tmp;
					eS[a][EII.init] += tmp * EII.ionB[i];

					mxw_tmp = v_t * Elecs.sigmaBEB(e_t, EII.ionB[i], EII.kin[i], EII.occ[i]);
					Sp[a][EII.init] += mxw_tmp;
					eSp[a][EII.init] += mxw_tmp * EII.ionB[i];
					X[a][EII.fin[i]] += *(map_p[a][EII.init] + adams_n) * (mxw_tmp * Elecs.Np[m] + tmp * Elecs.N[m]);
					W[a][EII.init] += v_t * Elecs.sigmaBEBw1(e_t, EII.ionB[i], EII.kin[i], EII.occ[i]);
				}		
			}

			tmp_dEdt = 0;
			tmp_dNdt = 0;
			tmp_dNpdt = 0;
			tmp_dEpdt = 0;
			for (int i = 0; i < tot_p_size; i++) {// Free up last element in p[i] and erase the first one.
				A[a][i] -= Elecs.N[m] * (S[a][i] + Elecs.N[m] * Tbr[a][i] ) + Elecs.Np[m] * Sp[a][i]; 
				tmp = A[a][i] * *(map_p[a][i] + adams_n) + X[a][i];
				*(map_dpdt[a][i] + adams_n) = tmp;
				tmp_dNpdt += *(map_p[a][i] + adams_n) * ( Pht[a][i] - Gesc * Elecs.Np[m] );
				tmp_dEpdt += *(map_p[a][i] + adams_n) * ( ePht[a][i] - Gesc * Elecs.Ep[m] - Elecs.Np[m] * (eSp[a][i] + W[a][i]) );
				tmp_dNdt += *(map_p[a][i] + adams_n) * ( Aug[a][i] + Elecs.N[m] * (S[a][i] - Elecs.N[m] * Tbr[a][i]) + Elecs.Np[m] * Sp[a][i] );
				tmp_dEdt += *(map_p[a][i] + adams_n) * ( eAug[a][i] - Elecs.N[m] * (eS[a][i] - Elecs.N[m] * eTbr[a][i]) + Elecs.Np[m] * W[a][i] );

				for (int j = 1; j < p_size; j++) {
					*(map_p[a][i] + j - 1) = *(map_p[a][i] + j);
					*(map_dpdt[a][i] + j - 1) = *(map_dpdt[a][i] + j);
				}
				*(map_p[a][i] + adams_n) = 0;
				*(map_dpdt[a][i] + adams_n) = 0;
			}
			Elecs.dNdt[m] += Store[a].nAtoms * tmp_dNdt;
			Elecs.dEdt[m] += Store[a].nAtoms * tmp_dEdt;
			Elecs.dNpdt[m] += Store[a].nAtoms * tmp_dNpdt;
			Elecs.dEpdt[m] += Store[a].nAtoms * tmp_dEpdt;

			for (int i = 0; i < tot_p_size; i++) {// Correct p.
				tmp = *(map_p[a][i] + adams_n - 2);
				for (int j = 0; j < adams_n; j++) {
					tmp += Moulton_5[j] * *(map_dpdt[a][i] + adams_n - j - 1) * dt[m - j];
				}
				*(map_p[a][i] + adams_n - 1) = tmp;
				if (tmp < 0 || tmp > 1) {
					printf("\n i = %d \n", i);
					return m;
				}

				A[a][i] = 0.;
				X[a][i] = 0.;
				Pht[a][i] = 0.;
				ePht[a][i] = 0.;
				S[a][i] = 0.;
				eS[a][i] = 0.;
				Tbr[a][i] = 0.;
				eTbr[a][i] = 0.;
				Sp[a][i] = 0.;
				eSp[a][i] = 0.;
				W[a][i] = 0.;
			}
		}
    
    // Photo-secondary collision.
    /*
    tmp = log(1.5*Temperature*sqrt(Temperature/Constant::Pi*Elecs.N[m]));
    tmp *= Elecs.Np[m]*Elecs.N[m]/v_t*64*pow(Constant::Pi, 1.5)*Elecs.BettaInt(e_t/Temperature);
    Elecs.dEdt[m] += tmp;
    Elecs.dEpdt[m] -= tmp;
    */
		// Correct plasma.
		Elecs.E[m] = Elecs.E[m-1];
		Elecs.N[m] = Elecs.N[m-1];
		Elecs.Ep[m] = Elecs.Ep[m-1];
		Elecs.Np[m] = Elecs.Np[m-1];
		for (int j = 0; j < adams_n; j++) {
			Elecs.E[m] += Moulton_5[j] * Elecs.dEdt[m - j] * dt[m - j];
			Elecs.N[m] += Moulton_5[j] * Elecs.dNdt[m - j] * dt[m - j];
			Elecs.Ep[m] += Moulton_5[j] * Elecs.dEpdt[m - j] * dt[m - j];
			Elecs.Np[m] += Moulton_5[j] * Elecs.dNpdt[m - j] * dt[m - j];
		}
		
		e_t = Elecs.Ep[m] / Elecs.Np[m];
		v_t = sqrt(2*e_t);
		Gesc = 1.*v_t / store.R;
		Temperature = 2*Elecs.E[m]/3/Elecs.N[m];
		Elecs.SetMaxwellPF(Temperature);

		// Correct plasma.
		Elecs.dNdt[m] = 0.;
		Elecs.dEdt[m] = 0.;
		Elecs.dNpdt[m] = 0.;
		Elecs.dEpdt[m] = 0.;

		for (int a = 0; a < Store.size(); a++) {
			tot_p_size = Store[a].num_conf;
			
			//correct dpdt
			for (auto& rate: Store[a].Photo) {
				tmp = rate.val * f[m];
				A[a][rate.from] -= tmp;
				Pht[a][rate.from] += tmp;
				ePht[a][rate.from] += tmp * rate.energy;
				X[a][rate.to] += tmp * *(map_p[a][rate.from] + adams_n - 1);
			}
			for (auto& rate: Store[a].Fluor) {
				A[a][rate.from] -= rate.val;
				X[a][rate.to] += rate.val * *(map_p[a][rate.from] + adams_n - 1);
			}
			for (auto& rate: Store[a].Auger) {
				A[a][rate.from] -= rate.val;
				X[a][rate.to] += rate.val * *(map_p[a][rate.from] + adams_n - 1);
			}
			for (auto& EII: Store[a].EIIparams) {
				for (int i = 0; i < EII.fin.size(); i++) {
					tmp = Pix4*Elecs.MaxwellEII(EII.ionB[i], EII.kin[i], EII.occ[i]);
					mxw_tmp = Pi3x8 * Elecs.MaxwellPF(EII.ionB[i]) * tmp;
					Tbr[a][EII.fin[i]] += mxw_tmp;
					eTbr[a][EII.fin[i]] += mxw_tmp * EII.ionB[i];
					X[a][EII.init] += mxw_tmp * *(map_p[a][EII.fin[i]] + adams_n - 1) * Elecs.N[m] * Elecs.N[m];
					S[a][EII.init] += tmp;
					eS[a][EII.init] += tmp * EII.ionB[i];

					mxw_tmp = v_t * Elecs.sigmaBEB(e_t, EII.ionB[i], EII.kin[i], EII.occ[i]);
					Sp[a][EII.init] += mxw_tmp;
					eSp[a][EII.init] += mxw_tmp * EII.ionB[i];
					X[a][EII.fin[i]] += *(map_p[a][EII.init] + adams_n - 1) * (mxw_tmp * Elecs.Np[m] + tmp * Elecs.N[m]);
					W[a][EII.init] += v_t * Elecs.sigmaBEBw1(e_t, EII.ionB[i], EII.kin[i], EII.occ[i]);
				}		
			}

			tmp_dNdt = 0;
			tmp_dEdt = 0;
			tmp_dNpdt = 0;
			tmp_dEpdt = 0;
			for (int i = 0; i < tot_p_size; i++) {
				A[a][i] -= Elecs.N[m] * (S[a][i] + Elecs.N[m] * Tbr[a][i] ) + Elecs.Np[m] * Sp[a][i]; 
				tmp = A[a][i] * *(map_p[a][i] + adams_n - 1) + X[a][i];
				*(map_dpdt[a][i] + adams_n - 1) = tmp;
				tmp_dNpdt += *(map_p[a][i] + adams_n - 1) * ( Pht[a][i] - Gesc * Elecs.Np[m] );
				tmp_dEpdt += *(map_p[a][i] + adams_n - 1) * ( ePht[a][i] - Gesc * Elecs.Ep[m] - Elecs.Np[m] * (eSp[a][i] + W[a][i]) );
				tmp_dNdt += *(map_p[a][i] + adams_n - 1) * ( Aug[a][i] + Elecs.N[m] * (S[a][i] - Elecs.N[m] * Tbr[a][i]) + Elecs.Np[m] * Sp[a][i] );
				tmp_dEdt += *(map_p[a][i] + adams_n - 1) * ( eAug[a][i] - Elecs.N[m] * (eS[a][i] - Elecs.N[m] * eTbr[a][i]) + Elecs.Np[m] * W[a][i] );
			}
			Elecs.dNdt[m] += Store[a].nAtoms * tmp_dNdt;
			Elecs.dEdt[m] += Store[a].nAtoms * tmp_dEdt;
			Elecs.dNpdt[m] += Store[a].nAtoms * tmp_dNpdt;
			Elecs.dEpdt[m] += Store[a].nAtoms * tmp_dEpdt;
		}
    
    // Photo-secondary collision.
    /*
    tmp = log(1.5*Temperature*sqrt(Temperature/Constant::Pi/Elecs.N[m]));
    tmp *= Elecs.Np[m]*Elecs.N[m]*64*pow(Constant::Pi, 1.5)*Elecs.BettaInt(e_t/Temperature)/v_t;
    Elecs.dEdt[m] += tmp;
    Elecs.dEpdt[m] -= tmp;
    */

		if ((m % storage_count) == 0) {
			int k = m / storage_count;
			time_storage[k] = t[m];
			for (int i = 0; i < p.size(); i++)	{
				p_storage[i][k] = p[i][p_size - 2];
			}
		}
	}

	cout << endl;
	return 0;
}

IntegrateRateEquation::~IntegrateRateEquation()
{
}

inline bool CompareChar(vector<char> &pat, char X)
{
	for (auto& v : pat)
	{
		if (v == X) return true;
	}
	return false;
}