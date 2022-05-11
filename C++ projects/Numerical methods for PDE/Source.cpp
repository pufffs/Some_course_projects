#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include<fstream>
#include <cassert>
#include <chrono>
using namespace std;
using namespace chrono;
class CB
{
public:
	CB(double T_, double F_, double R_, double r_, double kappa_, double mu_, double X_,
		double C_, double alpha_, double beta_, double sigma_, double Smax_, int imax_, int jmax_) :
		T(T_), F(F_), R(R_), r(r_), kappa(kappa_),
		mu(mu_), X(X_), C(C_), alpha(alpha_), beta(beta_), sigma(sigma_), Smax(Smax_), imax(imax_), jmax(jmax_){};
	//constructor.
	vector<double> EU_CN()
	{
		vector<double>vold(jmax + 1), vnew(jmax + 1);
		for (int j = 0; j <= jmax; j++)
		{
			vold[j] = max(F, R * j * dS);
		}
		for (int i = imax-1; i >= 0; i--)
		{
			double t = (i + 0.5) * dt; double t1 = i * dt;
			vector<vector<double>>M = matrix(vold, t, t1);
			vnew = thomasSolve(M[0], M[1], M[2], M[3]);
			vold = vnew;
		}
		return vnew;
	}

	vector<double> AMP(double P,double t0,double rho,double tol,int iterMax)
	{
		vector<double>vold(jmax + 1), vnew(jmax + 1);
		for (int j = 0; j <= jmax; j++)
		{
			vold[j] = max(F, R * j * dS);
		}
		for (int i = imax - 1; i >= 0; i--)
		{
			double t = (i + 0.5) * dt; double t1 = i * dt;
			if (t1 >= t0)
			{
				int penaltyIt;
				vector<vector<double>>M = matrix1(vold, t, t1);
				for (penaltyIt = 0; penaltyIt < iterMax; penaltyIt++)
				{
					vector<double> aHat(M[0]), bHat(M[1]), cHat(M[2]), dHat(M[3]);
					for (int j = 1; j < jmax; j++)
					{
						if (vnew[j] < R*dS*j)
						{
							bHat[j] = M[1][j] - rho; dHat[j] = M[3][j] - rho * (R * dS * j);
						}
					}
					vector<double> y = thomasSolve(aHat, bHat, cHat, dHat);
					double error = 0.;
					for (int j = 0; j <= jmax; j++)
						error += (vnew[j] - y[j]) * (vnew[j] - y[j]);
					vnew = y;
					if (error < tol * tol)
					{
						break;
					}
				}
				if (penaltyIt >= iterMax)
				{
					cout << " Error NOT converging within required iterations" << endl;
					throw;
				}
				vold = vnew;
			}
			if (t1 < t0)
			{
				{
					int penaltyIt;
					vector<vector<double>>M = matrix1(vold, t, t1);
					for (penaltyIt = 0; penaltyIt < iterMax; penaltyIt++)
					{
						vector<double> aHat(M[0]), bHat(M[1]), cHat(M[2]), dHat(M[3]);
						for (int j = 0; j < jmax; j++)
						{
							double D = max(R * dS * j, P);
							if (vnew[j] < D)
							{
								bHat[j] = M[1][j] - rho; dHat[j] = M[3][j] - rho * D;
							}

						}
						vector<double> y = thomasSolve(aHat, bHat, cHat, dHat);
						double error = 0.;
						for (int j = 0; j <= jmax; j++)
							error += (vnew[j] - y[j]) * (vnew[j] - y[j]);
						vnew = y;
						if (error < tol * tol)
						{
							break;
						}
					}
					if (penaltyIt >= iterMax)
					{
						cout << " Error NOT converging within required iterations" << endl;
						throw;
					}
					vold = vnew;
				}
			}
		}
		return vnew;
	}
	vector<double> AMPolicy(double P, double t0, double rho, double tol, int iterMax)
	{
		vector<double>vold(jmax + 1), vnew(jmax + 1);
		for (int j = 0; j <= jmax; j++)
		{
			vold[j] = max(F, R * j * dS);
		}
		for (int i = imax - 1; i >= 0; i--)
		{
			double t = (i + 0.5) * dt; double t1 = i * dt;
			if (t1 >= t0)
			{
				int penaltyIt;
				vector<vector<double>>M = matrix1(vold, t, t1);
				for (penaltyIt = 0; penaltyIt < iterMax; penaltyIt++)
				{
					vector<double> aHat(M[0]), bHat(M[1]), cHat(M[2]), dHat(M[3]);
					for (int j = 1; j < jmax; j++)
					{
						if (vnew[j] < R * dS * j)
						{
							aHat[j] = 0.; cHat[j] = 0.;
							bHat[j] = 1.;dHat[j] = R * dS * j;
						}
					}
					vector<double> y = thomasSolve(aHat, bHat, cHat, dHat);
					double error = 0.;
					for (int j = 0; j <= jmax; j++)
						error += (vnew[j] - y[j]) * (vnew[j] - y[j]);
					vnew = y;
					if (error < tol * tol)
					{
						break;
					}
				}
				if (penaltyIt >= iterMax)
				{
					cout << " Error NOT converging within required iterations" << endl;
					throw;
				}
				vold = vnew;
			}
			if (t1 < t0)
			{
				{
					int penaltyIt;
					vector<vector<double>>M = matrix1(vold, t, t1);
					for (penaltyIt = 0; penaltyIt < iterMax; penaltyIt++)
					{
						vector<double> aHat(M[0]), bHat(M[1]), cHat(M[2]), dHat(M[3]);
						for (int j = 0; j < jmax; j++)
						{
							double D = max(R * dS * j, P);
							if (vnew[j] < D)
							{
								aHat[j] = 0.; cHat[j] = 0.;
								bHat[j] = 1.; dHat[j] = D;
							}

						}
						vector<double> y = thomasSolve(aHat, bHat, cHat, dHat);
						double error = 0.;
						for (int j = 0; j <= jmax; j++)
							error += (vnew[j] - y[j]) * (vnew[j] - y[j]);
						vnew = y;
						if (error < tol * tol)
						{
							break;
						}
					}
					if (penaltyIt >= iterMax)
					{
						cout << " Error NOT converging within required iterations" << endl;
						throw;
					}
					vold = vnew;
				}
			}
		}
		return vnew;
	}
	double AMexplicit(double S0, double X, double T, double r, double sigma, int iMax, int jMax, double S_max)
	{
		double dS = S_max / jMax;
		double dt = T / iMax;
		vector<double> S(jMax + 1), vOld(jMax + 1), vNew(jMax + 1);
		for (int j = 0; j <= jMax; j++)
		{
			S[j] = j * dS;
		}
		for (int j = 0; j <= jMax; j++)
		{
			vOld[j] = max(X - S[j], 0.);
			vNew[j] = max(X - S[j], 0.);
		}
		for (int i = iMax - 1; i >= 0; i--)
		{
			vNew[0] = X * exp(-r * (T - i * dt));
			for (int j = 1; j <= jMax - 1; j++)
			{
				double A, B, C;
				A = 0.5 * sigma * sigma * j * j * dt + 0.5 * r * j * dt;
				B = 1. - sigma * sigma * j * j * dt;
				C = 0.5 * sigma * sigma * j * j * dt - 0.5 * r * j * dt;
				vNew[j] = 1. / (1. + r * dt) * (A * vOld[j + 1] + B * vOld[j] + C * vOld[j - 1]);
			}
			vNew[jMax] = 0.;
			vOld = vNew;
		}
		int jstar;
		jstar = S0 / dS;
		double sum = 0.;
		sum = sum + (S0 - S[jstar + 1]) / (S[jstar] - S[jstar + 1]) * vNew[jstar];
		sum = sum + (S0 - S[jstar]) / (S[jstar + 1] - S[jstar]) * vNew[jstar + 1];
		return sum;
	}
	double inter(double S0)
	{
		int js;
		js = S0 / dS;
		vector<double> v = EU_CN();
		vector<double> S1 = S();
		return  ( S0 - S1[js + 1] ) / ( S1[js] - S1[js + 1] ) * v[js] + ( S0 - S1[js] ) / ( S1[js + 1] - S1[js] ) * v[js + 1];
	}
	double Hinter(double S0)
	{
		int js;
		js = S0 / dS;
		vector<double> v = EU_CN();
		vector<double> S1 = S();
		double A = ((S0 - S1[js]) * (S0 - S1[js + 1])) / ((S1[js - 1] - S1[js]) * (S1[js - 1] - S1[js + 1])) * v[js - 1];
		double B = ((S0 - S1[js-1]) * (S0 - S1[js + 1])) / ((S1[js] - S1[js-1]) * (S1[js] - S1[js + 1])) * v[js];
		double C = ((S0 - S1[js-1]) * (S0 - S1[js])) / ((S1[js+1] - S1[js-1]) * (S1[js+1] - S1[js])) * v[js + 1];
		return A + B + C;
	}
	double inter1(double S0, double P, double t0, double rho, double tol, int iterMax)
	{
		int js;
		js = S0 / dS;
		vector<double> v = AMPolicy(P, t0, rho, tol, iterMax);
		vector<double> S1 = S();
		return  (S0 - S1[js + 1]) / (S1[js] - S1[js + 1]) * v[js] + (S0 - S1[js]) / (S1[js + 1] - S1[js]) * v[js + 1];
	}
	double Hinter1(double S0, double P, double t0, double rho, double tol, int iterMax,int c)
	{
		int js;
		js = S0 / dS;
		vector<double> v;
		if (c == 1)
		{
			v = AMP(P, t0, rho, tol, iterMax);
		}
		if (c == 2)
		{
			v = AMPolicy(P, t0, rho, tol, iterMax);
		}
		vector<double> S1 = S();
		double A = ((S0 - S1[js]) * (S0 - S1[js + 1])) / ((S1[js - 1] - S1[js]) * (S1[js - 1] - S1[js + 1])) * v[js - 1];
		double B = ((S0 - S1[js - 1]) * (S0 - S1[js + 1])) / ((S1[js] - S1[js - 1]) * (S1[js] - S1[js + 1])) * v[js];
		double C = ((S0 - S1[js - 1]) * (S0 - S1[js])) / ((S1[js + 1] - S1[js - 1]) * (S1[js + 1] - S1[js])) * v[js + 1];
		return A + B + C;
	}
	double check(double S1, double t1)
	{
		return R * S1 * exp(-(kappa + r) * (T - t1)) + X * R * exp(-r * (T - t1)) * (1. - exp(-kappa * (T - t1))) +
			C / (alpha + r) * (exp(-alpha * t1) - exp(-(alpha + r) * T) * exp(r * t1));
	}
private:
	double theta(double t)
	{
		return (1. + mu) * X * exp(mu * t);
	}
	double K(double t)
	{
		return C * exp(-alpha * t);
	}
	vector<double> S()
	{
		vector<double>v(jmax + 1);
		for (int j = 0; j <= jmax; j++)
		{
			v[j] = j * dS;
		}
		return v;
	}
	vector<vector<double>> matrix(const vector<double> &vold,double t,double t1)
	{
		vector<double> a(jmax + 1), b(jmax + 1), c(jmax + 1), d(jmax + 1);
		a[0] = 0.;
		b[0] = -1. / dt - kappa * theta(t1) / dS - r;
		c[0] = kappa * theta(t1) / dS;
		d[0] = -(1. / dt) * vold[0] - K(t1);
		for (int j = 1; j <= jmax - 1; j++)
		{
			a[j] = 0.25*(sigma * sigma * pow(j, 2. * beta) * pow(dS, 2. * (beta - 1.)) - kappa * theta(t) / dS + kappa * j);
			b[j] = -0.5 * sigma * sigma * pow(j, 2. * beta) * pow(dS, 2. * (beta - 1.)) - 1. / dt - r / 2.;
			c[j]= 0.25*(sigma * sigma * pow(j, 2. * beta) * pow(dS, 2. * (beta - 1.)) + kappa * theta(t) / dS - kappa * j);
			d[j] = -a[j] * vold[j - 1] - c[j] * vold[j + 1] - (b[j] + 2. / dt) * vold[j] - K(t);
		}
		a[jmax] = 0.;
		b[jmax] = 1.;
		c[jmax] = 0.;
		d[jmax] = R * Smax * exp(-(kappa + r)*(T - t1)) + X * R * exp(-r * (T - t1)) * ( 1. - exp(-kappa * (T - t1)) ) +
			C / (alpha + r) * ( exp(-alpha * t1) - exp(-(alpha + r) * T) * exp(r * t1) );
		return { a,b,c,d };
	}
	vector<vector<double>> matrix1(const vector<double>& vold, double t, double t1)
	{
		vector<double> a(jmax + 1), b(jmax + 1), c(jmax + 1), d(jmax + 1);
			a[0] = 0.;
			b[0] = -1. / dt - kappa * theta(t1) / dS - r;
			c[0] = kappa * theta(t1) / dS;
			d[0] = -(1. / dt) * vold[0] - K(t1);
			for (int j = 1; j <= jmax - 1; j++)
			{
				a[j] = 0.25 * (sigma * sigma * pow(j, 2. * beta) * pow(dS, 2. * (beta - 1.)) - kappa * theta(t) / dS + kappa * j);
				b[j] = -0.5 * sigma * sigma * pow(j, 2. * beta) * pow(dS, 2. * (beta - 1.)) - 1. / dt - r / 2.;
				c[j] = 0.25 * (sigma * sigma * pow(j, 2. * beta) * pow(dS, 2. * (beta - 1.)) + kappa * theta(t) / dS - kappa * j);
				d[j] = -a[j] * vold[j - 1] - c[j] * vold[j + 1] - (b[j] + 2. / dt) * vold[j] - K(t);
			}
			a[jmax] = 0.;
			b[jmax] = 1.;
			c[jmax] = 0.;
			d[jmax] = R * Smax;
			return { a,b,c,d };
	}
	vector<double> thomasSolve(const std::vector<double>& a, const std::vector<double>& b_, const std::vector<double>& c, std::vector<double>& d)
	{
		int n = a.size();
		std::vector<double> b(n), temp(n);
		// initial first value of b
		b[0] = b_[0];
		for (int j = 1; j < n; j++)
		{
			b[j] = b_[j] - c[j - 1] * a[j] / b[j - 1];
			d[j] = d[j] - d[j - 1] * a[j] / b[j - 1];
		}
		// calculate solution
		temp[n - 1] = d[n - 1] / b[n - 1];
		for (int j = n - 2; j >= 0; j--)
			temp[j] = (d[j] - c[j] * temp[j + 1]) / b[j];
		return temp;
	}

	double T, F, R, r, kappa, mu, X, C, alpha, beta, sigma, Smax;
	int imax, jmax;
	double dt = T / imax; double dS = Smax / jmax; double dt1 = T / (imax * imax);
};

int main()
{
	/*
	CB example1(5., 140., 2., 0.0202, 0.05, 0.0186, 67.58, 0.705, 0.01, 0.787, 0.625, 280., 500, 1000);
	CB example2(5., 140., 2., 0.0202, 0.05, 0.0186, 67.58, 1.41, 0.01, 0.787, 0.625, 280., 500, 1000);
	CB example3(5., 140., 2., 0.0202, 0.05, 0.0186, 67.58, 2.115, 0.01, 0.787, 0.625, 280., 500, 1000);
	vector<double> a1 = example1.AMP(150., 1.8849, 1e6, 1e-6, 1000);
	vector<double> a2 = example2.AMP(150., 1.8849, 1e6, 1e-6, 1000);
	vector<double> a3 = example3.AMP(150., 1.8849, 1e6, 1e-6, 1000);
	ofstream w("vs.csv");
	double ds = 280. / 1000.;
	for (int j = 0; j <= 1000; j++)
	{
		w << ds * j << ',' << a2[j]-a1[j] << ','<< a3[j]-a1[j] << endl;
		//w << ds * j << ',' << a1[j] << ','<<a2[j]<<',' << a3[j]<< endl;
	}
	*/
	/*
	CB example2(5., 140., 2., 0.0202, 0.05, 0.0186, 67.58, 1.41, 0.01, 0.787, 0.625,280., 4000, 1000);
	
	auto start = high_resolution_clock::now();
	double a=example2.Hinter1(67.58, 150., 1.8849, 1e6, 1e3, 1000,2);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout.precision(12); cout << "value :" << fabs(162.944801-a) << endl;
	cout << "Time taken:" << duration.count() << endl;
	*/
	/*
	ofstream w("out.csv");
	for (int n = 100; n <= 2000; n +=100)
	{
		CB example1(5., 140., 2., 0.0202, 0.05, 0.0186, 67.58, 1.41, 0.01, 0.787, 0.625, 280., n, 5000);
		CB example2(5., 140., 2., 0.0202, 0.05, 0.0186, 67.58, 1.41, 0.01, 0.787, 0.625, 280., 5000,n);
		double a1 = example1.inter1(67.58, 150., 1.8849, 1e6, 1e-6, 1000);
		double a2 = example2.inter1(67.58, 150., 1.8849, 1e6, 1e-6, 1000);
		w << n << ','<<fabs(162.944801 - a1) << ',' << fabs(162.944801 - a2) << endl;
	}
	*/
	
	double vo = 0.;
	double extravalue = 0.;
	auto start = high_resolution_clock::now();
	for (int n = 4; n <= 5000; n *= 2)
	{
		CB example2(5., 140., 2., 0.0202, 0.05, 0.0186, 67.58, 1.41, 0.01, 0.787, 0.625, 300.,n,500+n/20.);
		double vn = example2.inter1(67.58, 150., 1.8849, 1e6, 1e3, 1000);
		if (n >= 8)
		{
			extravalue = (4. * vn - vo) / 3.;
		}
		double dif = fabs(162.944801 - extravalue);
		if (dif < 1e-3)
		{
			cout.precision(12); cout << dif<< endl;
		}
		vo = vn;
	}
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	cout << "Time taken:" << duration.count() << endl;
	
	return 0;
}