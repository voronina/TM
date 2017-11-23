/*“естовый пример Ѕѕ‘:
procedure PPT(var x:mas;sign, NP:integer; T:real);
var
nmax,i,nn,mm,lr,nw,ii,j,loc,nw1:integer;
zz,W,delta:real;
MSK:array[1..20] of integer;
cs:array[1..2] of real;
cxcs;xa,hold:complex;
label 1,2;*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define L 64
#define N 4096
struct complex {
	float re, im;
};

float M_PI = acos(1.0);

void fft(struct complex x[], int sign, int np, float t);
struct complex cmplx(float x, float y);
struct complex add(struct complex x, struct complex y);
struct complex red(struct complex x, struct complex y);
struct complex mult(struct complex x, struct complex y);
int st(int n, int np);
float absc(struct complex x);

float test_function(float x) {
	if (x >= -8. && x < -4) {
		return(0.2);
	}
	else if (x >= -4 && x < 1) {
		return(1.0);
	}
	else if (x >= 1 && x <= 9) {
		return(-0.1);
	}
	else {
		return(0);
	}
}



int main(int argc, char* argv[]) {
/*
	int n[N];
	float xn[N];
	struct complex yn[N];
	float Wn[N];
	int j;
	float v = (float)L / (float)N;


	// ???
	for (j = 0; j < N; j++) {
		n[j] = j;
		xn[j] = (j - (N / 2))*v;
		yn[j].re = test_function(xn[j]);
		yn[j].im = 0;
		Wn[j] = (float)n[j] / (float)L;
	}
	// ???
	
	FILE * pFile;
	pFile = fopen("results.txt", "w");

	float remas[N];
	int i;
	fft(yn, 1, 12, (float)L);
	for (i = 0; i<N; i++) {
		remas[i] = sqrt(yn[i].re*yn[i].re + yn[i].im*yn[i].im);
		// printf("%f",remas[i]);
		fprintf(pFile, "%f %f\n", Wn[i], remas[i]);
	}
	fclose(pFile);
	*/
	int exponent = 6;
	int length = pow(2, exponent);
	float x_begin = -10;
	float x_end = 10;
	float step = (x_end - x_begin) / length;

	complex* series = (complex*) malloc(length*sizeof(complex));
	float* amplitudes = (float *) malloc(length*sizeof(float));
	int i;
	for (i = 0; i < length; i++) {
		series[i].re = test_function(x_begin + i*step);
		series[i].im = 0;
	}

	fft(series, 1, exponent, 1);

	/*for (i = 0; i < length / 2; i++) 
	{
		amplitudes[i] = absc(series[length / 2 + i]); 
		amplitudes[i + length / 2] = absc(series[i]); 
	}*/

	for (i = 0; i < length; i++)
	{
		amplitudes[i] = absc(series[i]); //???
	}

	FILE * pFile;
	pFile = fopen("results.txt", "w");

	for (i = 0; i < length; i++) {
		fprintf(pFile, "%d %f\n", i /* - length / 2*/, amplitudes[i]);
		//fprintf(pFile, "%d %f\n", i, series[i].re);
	}

	fclose(pFile);

	printf("e\n");
	//getchar();
}

float absc(struct complex x) {
	return sqrt(x.re*x.re + x.im*x.im);
}

void fft(struct complex x[], int sign, int np, float t) {
	int msk[19];
	float cs[2];
	struct complex cxcs, xa, hold;
	int nmax, i, j, nn, lr, nw, nw1, ii, ij, loc, ll, mm;
	float zz, pi, delta, w;
	pi = M_PI;
	nmax = st(2, np);
	printf("%d\n", nmax);
	
	zz = 2 * pi*sign / nmax;
	delta = t / nmax;
	if (sign<0) {
		delta = 1 / t;
	}
	msk[0] = nmax / 2;
	for (i = 1; i<np; i++) {
		msk[i] = msk[i - 1] / 2;
	}
	nn = nmax;
	mm = 2;
	for (lr = 1; lr <= np; lr++) {
		nn = nn / 2;
		nw = 0;
		for (i = 1; i <= mm; i = i + 2) {
			ii = nn*i;
			w = nw*zz;
			cs[0] = cos(w);
			cs[1] = sin(w);
			cxcs = cmplx(cs[0], cs[1]);
			for (j = 1; j <= nn; j++) {
				ii = ii + 1;
				ij = ii - nn;
				xa = mult(cxcs, x[ii - 1]);
				x[ii - 1] = red(x[ij - 1], xa);
				x[ij - 1] = add(x[ij - 1], xa);
			}
			for (loc = 2; loc <= np; loc++) {
				ll = nw - msk[loc - 1];
				if (ll <= 0) {
					break;
				}
				else {
					nw = ll;
				}
			}
			if (ll == 0) {
				nw = msk[loc];
			}
			else {
				nw = nw + msk[loc - 1];
			}
		}
		mm = 2 * mm;
	}
	nw = 0;
	for (i = 1; i <= nmax; i++) {
		nw1 = nw + 1;
		hold = x[nw1 - 1];
		if (nw1 - i>0) {
			x[nw1 - 1].re = x[i - 1].re*delta;
			x[nw1 - 1].im = x[i - 1].im*delta;
		}
		if (nw1 - i >= 0) {
			x[i - 1].re = hold.re*delta;
			x[i - 1].im = hold.im*delta;
		}
		for (loc = 1; loc <= np; loc++) {
			ll = nw - msk[loc - 1];
			if (ll <= 0) {
				break;
			}
			else {
				nw = ll;
			}
		}
		if (ll == 0) {
			nw = msk[loc];
		}
		else {
			nw = nw + msk[loc - 1];
		}
	}
}

struct complex cmplx(float x, float y) {
	struct complex z;
	z.re = x;
	z.im = y;
	return(z);
}

struct complex add(struct complex x, struct complex y) {
	struct complex z;
	z.re = x.re + y.re;
	z.im = x.im + y.im;
	return(z);
}

struct complex red(struct complex x, struct complex y) {
	struct complex z;
	z.re = x.re - y.re;
	z.im = x.im - y.im;
	return(z);
}

struct complex mult(struct complex x, struct complex y) {
	struct complex z;
	z.re = x.re*y.re - x.im*y.im;
	z.im = x.re*y.im + x.im*y.re;
	return(z);
}

int st(int n, int np) {
	int m, i;
	m = 1;
	for (i = 1; i <= np; i++) {
		m = m*n;
	}
	return(m);
}
