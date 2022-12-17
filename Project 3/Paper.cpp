#define _USE_MATH_DEFINES
#include<iomanip>
#include<iostream>
#include<math.h>
#include<cmath>
using namespace std;
class Alpha {
private:
	int Z1, l, N, Q;
	long double C1, C2,M,y,r1,r2,R;
	const long double e = 1.439967;
    const double u = 931.5;
	const double h = 197.326980;
	long double phi1,phi2;
	long double r = 0.00000000001;
	long double p;

public:
	Alpha(int x1, int x2, int x3, int x4) {
		Z1 = x1;
		N = x2;
		l = x3;
		Q = x4;
		long double m1 = (Z1 + N - 4);long double m2 = 4;
		C1 = 1.28 * pow(m1, 1.0 / 3) - 0.76 + 0.8 * pow(m1, -1.0 / 3) - pow(1.28 * pow(m1, 1.0 / 3) - 0.76 + 0.8 * pow(m1, -1.0 / 3), -1);
	    C2 = 1.28 * pow(m2, 1.0 / 3) - 0.76 + 0.8 * pow(m2, -1.0 / 3) - pow(1.28 * pow(m2, 1.0 / 3) - 0.76 + 0.8 * pow(m2, -1.0 / 3), -1);
        M = ((Z1 + N - 4.0) * u * 4.0 * u) / ((Z1 + N - 4) * u + 4.0 * u);
	    y = 0.9517 * (1 - 1.7826 * pow((N - Z1), 2) / (double)pow((N + Z1), 2));
	}
	long double Vc() {
		if (r < r1)
			Z1 * 2 * pow(e, 2) * pow(2 * R, -1) * (3 - pow((double)r / R, 2));
		else
		return ((Z1-2) *2 * e) / (double)r;
	}
	long double Va(){
		return (pow(h, 2) * l * (l + 0.5)) / double(2 * M * pow(r, 2));
		
	}
	long double Vp() {
		long double z = r - C1 - C2;
		if (z > 1.9475) {
			phi1 = -4.41 * pow(M_E, z / (-0.7176));
			long double c= 4 * M_PI * y * ((C1 * C2) / (double)(C1 + C2)) * phi1;
			return c;
		}
		else if (z > 0 && z <= 1.9745) {
			phi2 = -1.7817 + 0.927 * z + 0.0169 * pow(z, 2) + -0.05148 * pow(z, 3);
			long double c= 4 * M_PI * y * ((C1 * C2) / (C1 + C2)) * phi2;
			return c;
		}
	}
	long double V(long double x) {
		r = x;
		return Va()+Vp()+Vc()-Q;
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
void root(Alpha I) {
	double x0, x1, x, f0, f1, f, e;
	x0 = 1;
	int step = 1;
	x1 = 100;
	e = 5 * pow(10, -4);
	do {
		f0 = I.V(x0++);
		f1 = I.V(x1);
		cout << f0 << endl;
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
		cout << x;
	}
		int main() {
			Alpha I(55, 57, 2, 0.823);
			root(I);

		}
