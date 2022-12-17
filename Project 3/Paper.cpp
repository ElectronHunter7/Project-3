#define _USE_MATH_DEFINES
#include<iomanip>
#include<iostream>
#include<math.h>
#include<cmath>
using namespace std;
class Alpha {
private:
	int Z1, l, N, Q;
	long double C1, C2, c1, c2, M, y, r1, r2, R, I1, I2,t1;
	const long double e = 1.439967;
    const double u = 931.5;
	const double h = 197.326980;
	long double phi1,phi2;
	long double r = 0.1;
	long double p;

public:
	Alpha(int x1, int x2, int x3, int x4) {
		Z1 = x1;
		N = x2;
		l = x3;
		Q = x4;
		int A1 = (Z1 + N - 1);int A2 = 1;
		I1 = (N - Z1 - 1) * pow(A1, -1);
		I2 = -1;
		r1 = 1.2332 * pow(A1, 1.0 / 3) + 2.8961 * pow(pow(A1, 2.0 / 3), -1) - 0.18688 * pow(A1, 1.0 / 3) * I1;
		r2 = 1.2332 * pow(A2, 1.0 / 3) + 2.8961 * pow(pow(A2, 2.0 / 3), -1) - 0.18688 * pow(A2, 1.0 / 3) * I2;
		c1 = r1 * (1 - 7 * pow(2 * pow(r1, 2), -1) - 49 * pow(8 * pow(r1, 4), -1));
		c2 = r2 * (1 - 7 * pow(2 * pow(r2, 2), -1) - 49 * pow(8 * pow(r2, 4), -1));
		t1 = 3 / 2 * (1.14) * ((32.65) * I1 - (1.0 / 12) * 0.75789 * (Z1 - 1) * pow(A1, -1.0 / 3) * pow(Q +9.0 / 4 + 32.65 * pow(A1, -1.0 / 3), -1));
		C1 = c1 + (N / A1) * t1;
		C2 = c2;
		R = (C1* C2)/(C1+C2);
		M = ((Z1 + N - 1) * u * 1 * u) / ((Z1 + N - 1) * u + 1 * u);
		y = 1.25284 * (1 - 2.345 * pow((N - Z1) / (N + Z1), 2));
	}
	long double Vc() {
		if (r < r1+r2)
		return (Z1-1) * e * pow(2 * (r1+r2), -1) * (3 - pow((double)r / (r1+r2), 2));
		else
		return ((Z1-1) * e) / (double)r;
	}
	long double Va(){
		return (pow(h, 2) * pow((l +0.5),2)) / double(2 * M * pow(r, 2));
		
	}
	long double Vp() {
		long double z = r - C1 - C2;
		if (z < 0) {
			phi1 = -1.7817 + 0.972 * z + 0.143 *pow(z,2) - 0.09 * pow(z, 3);
			return 4 * M_PI * y * R * phi1;
		}
		else if (0 < z && z < 1.9475) {
			phi1= -1.7817 + 0.9720 * z + 0.01696 * z * z - 0.05148 * z * z * z;
			return 4 * M_PI * y * R * phi1;
		}
		else if (z > 1.9475) {
			phi1 = -4.41 * exp(-z / (0.7176));
			return 4 * M_PI * y * R * phi1;
		}
	}
	long double V(long x) {
		r = x;
		return Va() + Vc() + Va() - Q;
	}
	void T() {
		long double c = pow(2 * M * (Vp() + Va() + Vc() - Q), 0.5);
		r = r1;
		while (r < r2) {
			r += 0.0001;
			c += c * r;
		}
		p = pow(M_E, (-2 / h) * c);
		long double v = 41 / (h * pow(Z1 + N, 1.0 / 3));
		cout<<log(2.0)/(v*p);
	}
	

		
};
double root(Alpha I) {
	double x0, x1, x, f0, f1, f, e;
	x0 = 1;
	int step = 1;
	x1 = 100;
	e = 5 * pow(10, -4);
	do {
		f0 = I.V(x0++);
		f1 = I.V(x1);
	}while (f0 * f1 > 0);
		do {
			x = x0 - (x0 - x1) * f0 * pow(f0 - f1, -1);
			f = I.V(x);
			if (f0 * f < 0) {
				x1 = x;
				f1 = f;

			}
			else {
				x0 = x;
				f0 = f;
			}
			step = step + 1;
		} while (fabs(f) > e);

		return x;
	}
int main() {
	Alpha I(69, 76, 5, 1.741);
	cout << root(I);

}
