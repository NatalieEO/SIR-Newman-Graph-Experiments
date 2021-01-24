
/* Author: Maleq Khan */


#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <cmath>
#include <ctime>
#include <cstdlib>

		
using namespace std;


// -------- initialize random number generator --------
void InitRand()
{
    srand(time(NULL));
}

void InitRand(unsigned int seed)
{
    srand(seed);
}


//  Return uniformly distributed double in range [0.0, 1.0].
inline double Uniform()
{
	return ((double)rand()) / (double) RAND_MAX;
}


//  Return uniformly distributed double in range [min, max).
double Uniform(double min, double max)
{
	return min + (max - min) * Uniform();
}


//  Return uniformly distributed integer in range [min, max].
int Uniform(int min, int max)
{
	return min + (int)floor((max - min + 1) * Uniform());
}

//  Return uniformly distributed integer in range [min, max].
unsigned long int Uniform(unsigned long int min, unsigned long int max)
{
	return min + (unsigned long int)floor((max - min + 1) * Uniform());
}


//  Return uniformly distributed integer in range [min, max].
char Uniform(char min, char max)
{
	return min + (char)floor((max - min + 1) * Uniform());
}


//  Return uniformly distributed integer in range [min, max].
unsigned int Uniform(unsigned int min, unsigned int max)
{
	return min + (unsigned int)floor((max - min + 1) * Uniform());
}


//  Return uniformly distributed integer in range [min, max].
unsigned char Uniform(unsigned char min, unsigned char max)
{
	return min + (unsigned char)floor((max - min + 1) * Uniform());
}



// Return a binomial random number with parameter n and p
// using inverse transform method -- efficient only if np is small.
int BinomialInv(int n, double p)
{
	double Pi, F, U, fp;
	int i;
	
	Pi = pow(1-p, n);   
	F = Pi;
	fp = p / (1-p);
	U = Uniform();
	for (i=0; F < U; i++) {
		Pi *= (fp * (n-i)) / (i+1);
		F += Pi;
	}
	return i;
}


// Return a binomial random number with parameter n and p
// using BPTE algorithm -- valid only if n*min(p, 1-p) > 10.
int BinomialBPTE(int n, double p)
{
    int y, M, k, i, x1, f1, z, w, x2, f2, z2, w2;
    double r, q, fM, p1, xM, xL, xR, c, a, lL, lR, p2, p3, p4;
    double u, v, x, s, F, rho, t, A;
        
	//Step 0.
    r = (p < 0.5) ? p : (1.0 - p);
    q = 1.0 - r;
    fM = n * r + r;
    M = (int) fM;
    p1 = ((int)(2.195*sqrt(n*r*q) - 4.6*q)) + 0.5;
    xM = M + 0.5;
    xL = xM - p1;
    xR = xM + p1;
    c = 0.134 + 20.5 / (15.3 + M);
    a = (fM - xL) / (fM - xL*r);
    lL = a * (1 + a / 2);
    a = (xR - fM) / (xR * q);
    lR = a * (1 + a / 2);
    p2 = p1 * (1 + 2*c);
    p3 = p2 + c / lL;
    p4 = p3 + c / lR;
	
	while (true) {
		// Step 1.
		u = Uniform(0.0, p4);
		v = Uniform();
		if (u <= p1) {
                    y = int(xM - p1*v + u);
                    break;
		}
		
		// Step 2.
		if (u <= p2) {
			x = xL + (u - p1) / c;
			v = v * c + 1 - abs(M-x+0.5) / p1;
			if (v > 1) continue;
			y = (int) x;
		}
		// Step 3.
		else if (u <= p3) {
			y = (int) (xL + log(v) / lL);
			if (y < 0) continue;
			v = v * (u - p2) * lL;
		}
		//Step 4.
		else {
			y = (int) (xR - log(v) / lR);
			if (y > n) continue;
			v = v * (u - p3) * lR;
		}
		// Step 5.
		k = abs(y-M);		
		if (k <= 20 || k >= (n*r*q)/2 - 1) {		
		//step 5.1		
			s = r / q;
			a = s * (n+1);
			F = 1.0;
			if (M < y) {
				i = M;
				do {
					i++;
					F *= a / i - s;
				} while (i < y);
			}
			else if (M > y) {
				i = y;
				do {
					i++;
					F /= a / i - s;
				} while (i < M);
			}
			if (v > F) continue;		// goto Step 1
			break;						// Otherwise, goto Step 6 
		}
		//step 5.2 -- Squeezing. Check the value of log(v) against the upper and
		// lower bound of log f(y)
		rho = (k/(n*r*q)) * ((k*(k/3.0 + 0.625) + 1.0/6.0) / (n*r*q) + 0.5);
		t = -(k*k) / (2*n*r*q);
		A = log(v);
		if (A < t - rho) break;			// goto Step 6
		if (A > t + rho) continue;		// goto Step 1
		
		//Step 5.3 -- final acceptance / rejection test
		x1 = y + 1;
		f1 = M + 1;
		z = n + 1 - M;
		w = n - y + 1;
		x2 = x1 * x1;
		f2 = f1 * f1;
		z2 = z * z;
		w2 = w * w;
		if (A > xM * log((double)f1/x1) 
				+ (n - M + 0.5)*log((double)z/w) 
				+ (y-M)*log((double)(w*r)/(x1*q)) 
				+ (13860 - (462 - (132 - (99 - 140.0/f2)/f2)/f2)/f2)/f1/166320
				+ (13860 - (462 - (132 - (99 - 140.0/z2)/z2)/z2)/z2)/z/166320
				+ (13860 - (462 - (132 - (99 - 140.0/x2)/x2)/x2)/x2)/x1/166320
				+ (13860 - (462 - (132 - (99 - 140.0/w2)/w2)/w2)/w2)/w/166320
		   ) continue;
	} // end while	

	//Step 6.
	if (p > 0.5) return (n - y);		
	return y;
        return 0;
}


// Return a binomial random number with parameter n and p
// using a combination of Inverse transform (Inv) and BPTE algorithm 
int Binomial(int n, double p)
{
    if (p == 0.0) return 0;
    if (p == 1.0) return n;
    if (n*p < 30)
	return BinomialInv(n, p);
    if (n*(1-p) < 30)
        return n - BinomialInv(n, 1-p);
    return BinomialBPTE(n, p);
}


// Generate k multinomial random variables N (a vector of k integers)
// using parameters n and parobability vector P of size k.
int *Multinomial(int n, double *P, int k, int *N)
{
	int i;
	double cp = 0.0;		// cumulative prob. - set in the next line

        for (i=0; i<k; i++)         // normalize vector P if the 
            cp += P[i];             // probabilties do not sum to 1.
	for (i=0; i<k-1; i++) {
            if (cp > 0) {
                N[i] = Binomial(n, P[i]/cp);
                n = n - N[i];
                cp = cp - P[i];
            }
            else N[i] = 0;
	}
	N[k-1] = n;                 // the last event of multinomial
	return N;
}



#endif // RANDOM_HPP
