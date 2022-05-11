#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <algorithm>
#include<iomanip>
#include<fstream>
#include <cassert>
using namespace std;
class Payoff
{
public:
	virtual double operator()(double ST) = 0;
	virtual double operator()(double ST, double A) = 0;
};
class EUROC :public Payoff
{
public:
	EUROC(double X_) : X(X_) {};
	virtual double operator()(double ST)
	{
		return max(ST - X, 0.);
	}
	virtual double operator()(double ST, double A) { return 0; }
private:
	double X;
};


class pricing
{
public:
	pricing(double T_, double sigma_, double r_, double S0_,double t_,double X_) :
		T(T_), sigma(sigma_), r(r_),S0(S0_),t(t_),X(X_) {}; //constructor
	double MC(Payoff&p,int N)
	{
		static mt19937 rnd;
		normal_distribution<> ND(0.0, 1.0);
		double sum = 0;
		for (int i = 0; i < N; i++)
		{
			double Z = ND(rnd);
			double ST= S0 * exp( (T - t) * (r -0.5 * sigma * sigma) + sigma * sqrt(T - t) * Z );
			sum += p(ST);
		}
		return exp(-r * (T - t)) * sum / N;
	}

	double AMC(Payoff& p, int N)
	{
		static mt19937 rnd;
		normal_distribution<> ND(0.0, 1.0);
		double sum = 0;
		for (int i = 0; i < N; i++)
		{
			double Z = ND(rnd);
			double ST_p = S0 * exp((T - t) * (r -0.5 * sigma * sigma) + sigma * sqrt(T - t) * Z);
			sum += p(ST_p);
			double ST_m= S0 * exp((T - t) * (r -0.5 * sigma * sigma) - sigma * sqrt(T - t) * Z);
			sum += p(ST_m);
		}
		return (exp(-r * (T - t)) * sum )/ (2. * N);
	}

	double MM(Payoff& p,int N)
	{
		static mt19937 ran;
		normal_distribution<> ND(0.0, 1.0);
		double sum = 0;
		vector<double> path;
		for (int i = 0; i < N; i++)
		{
			double Z = ND(ran);
			path.push_back(Z);
			sum += Z * Z;
		}
		double sd = sqrt(sum / N);
		double sum1 = 0;
		for (int i = 0; i < N; i++)
		{
			double phi = path[i]/sd;
			double ST_p = S0 * exp((T - t) * (r -0.5 * sigma * sigma) + sigma * sqrt(T - t) * phi);
			sum1 += p(ST_p);
			double ST_m = S0 * exp((T - t) * (r - 0.5 * sigma * sigma) - sigma * sqrt(T - t) * phi);
			sum1 += p(ST_m);
		}
		return (exp(-r * (T - t)) * sum1) / (2. * N);
	}
	double BO(double B, int K, int N)
	{
		double sum = 0.;
		for (int i = 0; i < N; i++)
		{
			vector<double> out = RW(B,K);
			sum += max(out[0]-X,0.);
		}
		return (exp(-r * (T - t)) * sum) / N;
	}

	double ABO(double B, int K, int N)
	{
		double sum = 0.;
		for (int i = 0; i < N; i++)
		{
			vector<double> out = ARW(K);
			if (out[1] <= B)
			{
				sum += max(out[0] - X, 0.);
			}
			if (out[3] <= B)
			{
				sum += max(out[2] - X, 0.);
			}
		}
		return (exp(-r * (T - t)) * sum) / (2.*N);
	}
	double MBO(double B, int K, int N)
	{
		double sum = 0.;
		for (int i = 0; i < N; i++)
		{
			vector<double> out = MRW(K);
			if (out[1] <= B)
			{
				sum += max(out[0] - X, 0.);
			}
			if (out[3] <= B)
			{
				sum += max(out[2] - X, 0.);
			}
		}
		return (exp(-r * (T - t)) * sum) / (2. * N);
	}


	double WHBO(double B, int m, int N)
	{
		double sum = 0.;
		for (int i = 0; i < N; i++) //loop over N simulations
		{
			vector<double> pair = WH(m); //exstract the pair (X,X^) 
			double current_position = pair[0];
			double running_maximum = pair[1];
			if (S0 * exp(T* (r- 0.5 * sigma * sigma) + sigma * running_maximum) <= B) 
			//Barrier condition check
			{
				double vs = S0 * exp(T * (r- 0.5 * sigma * sigma) 
				+ sigma * current_position)-X;   //payoff of call
				sum += max(vs, 0.);
			}
		}
		return (exp(-r *T) * sum) / N; //return option value
	}

	double Var(Payoff& p, int N, int M, int method)
	{
		double s = 0;
		double var = 0;
		vector<double> v(M);
		for (int m = 0; m < M; m++)
		{
			if (method == 0) { v[m] = MC(p, N); }
			if (method == 1) { v[m] = AMC(p, N); }
			if (method == 2) { v[m] = MM(p, N); }
			s += v[m];
		}
		double mean = s / M;
		for (int i = 0; i < M; i++)
		{
			var += pow(v[i] - mean, 2);
		}
		double svar = var / (M - 1.);
		return svar;
	}

	void CI(Payoff&p,int N,int M,int method)
	{
		double s = 0;
		double var = 0;
		vector<double> v(M);
		for (int m = 0; m < M; m++)
		{
			if (method==0){ v[m] = MC(p, N);}
			if (method==1){ v[m] = AMC(p, N);}
			if (method == 2) { v[m] = MM(p, N); }
			s += v[m];
		}
		double mean = s / M;
		for (int i = 0; i < M; i++)
		{
			var += pow(v[i] - mean, 2);
		}
		double svar = var / (M - 1.);
		//cout << "mean:" << mean << endl;
		//cout << "variance:" << svar << endl;
		double sd = sqrt(svar / M);
		cout << "$[" << mean - 2. * sd << "," << mean + 2. * sd << "]$";
	}

	vector<double> GBM(int N)
	{
		static mt19937 rnd;
		normal_distribution<> ND(0.0, 1.0);
		vector<double> record(N + 1);
		record[0] = S0;
		for (int i = 1; i <= N; i++)
		{
			double Z = ND(rnd);
			record[i] = record[i-1] * exp((T - t) * (r - 0.5 * sigma * sigma) + sigma * sqrt(T - t) * Z);
		}
		return record;
	}

	vector<double> BM(int N)
	{
		static mt19937 rnd;
		normal_distribution<> ND(0.0, 1.0);
		vector<double> record(N + 1);
		double dt = T / N;
		record[0] = 0.;
		for (int i = 1; i <= N; i++)
		{
			double Z = ND(rnd);
			record[i] = record[i - 1] * sqrt(dt) * Z;
		}
		return record;
	}

private:
	vector<double> WH(int n)
	{
		static mt19937 rnd;
		double lambda = sqrt(2. * n / T);
		exponential_distribution<double> Ex(lambda);
		vector<double> X(n + 1), Xb(n + 1);
		X[0] = Xb[0] = 0.;
		for (int i = 1; i <= n; i++)
		{
			double Sup = Ex(rnd);
			double Inf = -Ex(rnd);
			X[i] = X[i - 1] + Sup + Inf;
			Xb[i] = max(Xb[i - 1], X[i - 1] + Sup);
		}
		return { X[n],Xb[n] };
	}
	vector<double> RW(double B, int n)
	{
		static mt19937 rnd;
		normal_distribution<double> N(0, 1.);
		double dt = T / n;
		double a = sqrt(dt);
		vector<double> X(n + 1);
		X[0] = S0;
		double maxi = 0.;
		for (int i = 1; i <= n; i++)
		{
			const double Z = N(rnd);
			X[i] = X[i - 1] * exp(dt * (r - 0.5 * sigma * sigma) + sigma * a * Z);
			if (X[i] > maxi)
			{
				maxi = X[i];
			}
			if (maxi > B)
			{
				return { 0.,0. };
			}
		}
		return { X[n],maxi };
	}
	vector<double> ARW(int n)
	{
		static mt19937 rnd;
		normal_distribution<double> N(0, 1.);
		double dt = T / n;
		double a = sqrt(dt);
		vector<double> X(n + 1), X1(n + 1);
		X[0] = S0; X1[0] = S0;
		double maxi = 0.;
		double maxii = 0.;
		for (int i = 1; i <= n; i++)
		{
			const double Z = N(rnd);
			X[i] = X[i - 1] * exp(dt * (r - 0.5 * sigma * sigma) + sigma * a * Z);
			X1[i] = X1[i - 1] * exp(dt * (r - 0.5 * sigma * sigma) - sigma * a * Z);
			if (X[i] > maxi)
			{
				maxi = X[i];
			}
			if (X1[i] > maxii)
			{
				maxii = X1[i];
			}
		}
		return { X[n],maxi,X1[n],maxii };
	}
	vector<double> MRW(int n)
	{
		static mt19937 rnd;
		normal_distribution<double> N(0, 1.);
		double dt = T / n;
		double a = sqrt(dt);
		vector<double> path;
		vector<double> X(n + 1),X1(n + 1);
		X[0] = S0; X1[0] = S0;
		double maxi = 0.; 
		double maxii = 0.;
		double ss = 0.;
		for (int i = 0; i < n; i++)
		{
			const double Z = N(rnd);
			path.push_back(Z);
			ss += Z * Z;
		}
		double sd = sqrt(ss / n);
		for (int i = 1; i <= n; i++)
		{
			double phi = path[i - 1] / sd;
			X[i] = X[i - 1] * exp(dt * (r - 0.5 * sigma * sigma) + sigma * a * phi);
			X1[i] = X1[i - 1] * exp(dt * (r - 0.5 * sigma * sigma) - sigma * a * phi);
			if (X[i] > maxi)
			{
				maxi = X[i];
			}
			if (X1[i] > maxii)
			{
				maxii = X1[i];
			}
		}
		return { X[n],maxi,X1[n],maxii};
	}
	double T, sigma, r,S0, t;
	double X;
};
vector<double> WH1(double T, int n)
{
	static mt19937 rnd;
	double lambda = sqrt(2. * n / T);
	exponential_distribution<double> Ex(lambda);
	vector<double> X(n + 1), Xb(n + 1);
	X[0] = Xb[0] = 0.;
	for (int i = 1; i <= n; i++)
	{
		double Sup = Ex(rnd);
		double Inf = -Ex(rnd);
		X[i] = X[i - 1] + Sup + Inf;
		Xb[i] = max(Xb[i - 1], X[i - 1] + Sup);
	}
	return X;
}
vector<double> BM1(double T,int N)
{
	static mt19937 rnd;
	normal_distribution<> ND(0.0, 1.0);
	vector<double> record(N + 1);
	double dt = T / N;
	record[0] = 0.;
	for (int i = 1; i <= N; i++)
	{
		double Z = ND(rnd);
		record[i] = record[i - 1] * sqrt(dt) * Z;
	}
	return record;
}
vector<double> WH2(double T, int n)
{
	static mt19937 rnd;
	double lambda = sqrt(2. * n / T);
	exponential_distribution<double> Ex(lambda);
	vector<double> X(n + 1), Xb(n + 1);
	X[0] = Xb[0] = 0.;
	for (int i = 1; i <= n; i++)
	{
		double Sup = Ex(rnd);
		double Inf = -Ex(rnd);
		X[i] = X[i - 1] + Sup + Inf;
		Xb[i] = max(Xb[i - 1], X[i - 1] + Sup);
	}
	return {X[n],Xb[n]};
}
vector<double> BM2(double T, int N)
{
	static mt19937 rnd;
	normal_distribution<> ND(0.0, 1.0);
	vector<double> record(N + 1);
	double dt = T / N;
	record[0] = 0.;
	double maximum = 0.;
	for (int i = 1; i <= N; i++)
	{
		double Z = ND(rnd);
		record[i] = record[i - 1] * sqrt(dt) * Z;
		if (record[i] > maximum)
		{
			maximum = record[i];
		}
	}
	return {record[N],maximum};
}
int main()
{
	pricing a(1., 0.5, 0.05, 100., 0., 120.);
	double a2 = a.BO(220., 50, 5000);
	double a3 = a.BO(220., 50, 10000);
	double a4 = a.BO(220., 50, 50000);
	double a5 = a.BO(220., 50, 100000);
	cout.precision(5); cout<< (a2 - 5.88) / 5.88 * 100 <<"%" <<(a3 - 5.88) / 5.88 * 100 
		<<"%"  <<  (a4 - 5.88) / 5.88 * 100 << "%"  << (a5 - 5.88) / 5.88 * 100 <<"%"<<endl;
	
	double a6 = a.WHBO(220., 50, 5000);
	double a7 = a.WHBO(220., 50, 10000);
	double a8 = a.WHBO(220., 50, 50000);
	double a9 = a.WHBO(220., 50, 100000);
	cout.precision(5); cout << (a6 - 5.88) / 5.88 * 100 << "%" << (a7 - 5.88) / 5.88 * 100
		<< "%" << (a8 - 5.88) / 5.88 * 100 << "%" << (a9 - 5.88) / 5.88 * 100 << "%" << endl;
	return 0;
}