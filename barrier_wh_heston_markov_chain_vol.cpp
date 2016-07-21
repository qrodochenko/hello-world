#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include "pnl_fft.h"

#define SQR(X) (X*X)
#define MAX(A,B) ( (A) > (B) ? (A):(B))
#define MIN(A,B) ( (A) < (B) ? (A):(B))
#define MEMORY_ALLOCATION_FAILURE 0
#define OK 1
#define uint unsigned int 
/*Vectors and matrices*/

static double **V, **P_old, **P_new;
static double **y, **f;
static int **f_down, **f_up;
static int **y_down, **y_up;
static double **pu_y, **pd_y;
static double **pu_f, **pd_f;
static double ***S, ***F;
static double *ba_log_prices; /*basic asset price line, of length M*/
static double * fftfreqs; /*fft frequencies*/
/*Memory allocation*/
static int memory_allocation(uint Nt, uint N, uint M)
{
	uint i;

	/*V is the Nt+1 x Nt+1 matrice, storing volatility values*/
	V = (double**)calloc(Nt + 1, sizeof(double*));
	if (V == NULL)
		return MEMORY_ALLOCATION_FAILURE;
	for (i = 0; i<Nt + 1; i++)
	{
		V[i] = (double *)calloc(Nt + 1, sizeof(double));
		if (V[i] == NULL)
			return MEMORY_ALLOCATION_FAILURE;
	}
	/*pu_f is the Nt+1 x Nt+1 matrice*/
	pu_f = (double**)calloc(Nt + 1, sizeof(double*));
	if (pu_f == NULL)
		return MEMORY_ALLOCATION_FAILURE;
	for (i = 0; i<Nt + 1; i++)
	{
		pu_f[i] = (double *)calloc(Nt + 1, sizeof(double));
		if (pu_f[i] == NULL)
			return MEMORY_ALLOCATION_FAILURE;
	}
	/*pd_f is the Nt+1 x Nt+1 matrice*/
	pd_f = (double**)calloc(Nt + 1, sizeof(double*));
	if (pd_f == NULL)
		return MEMORY_ALLOCATION_FAILURE;
	for (i = 0; i<Nt + 1; i++)
	{
		pd_f[i] = (double *)calloc(Nt + 1, sizeof(double));
		if (pd_f[i] == NULL)
			return MEMORY_ALLOCATION_FAILURE;
	}

	/*f_down is the Nt+1 x Nt+1 matrice*/
	f_down = (int**)calloc(Nt + 1, sizeof(int*));
	if (f_down == NULL)
		return MEMORY_ALLOCATION_FAILURE;
	for (i = 0; i<Nt + 1; i++)
	{
		f_down[i] = (int *)calloc(Nt + 1, sizeof(int));
		if (f_down[i] == NULL)
			return MEMORY_ALLOCATION_FAILURE;
	}
	/*f_up is the Nt+1 x Nt+1 matrice*/
	f_up = (int**)calloc(Nt + 1, sizeof(int*));
	if (f_up == NULL)
		return MEMORY_ALLOCATION_FAILURE;
	for (i = 0; i<Nt + 1; i++)
	{
		f_up[i] = (int *)calloc(Nt + 1, sizeof(int));
		if (f_up[i] == NULL)
			return MEMORY_ALLOCATION_FAILURE;
	}
	/*P_old is the N+1 x N+1 matrice*/
	P_old = (double **)malloc((N + 1)*sizeof(double*));
	for (i = 0; i <= N; i++)
		P_old[i] = (double *)malloc((Nt + 1)*sizeof(double));
	/*P_new is the N+1 x N+1 matrice*/
	P_new = (double **)malloc((N + 1)*sizeof(double*));
	for (i = 0; i <= N; i++)
		P_new[i] = (double *)malloc((Nt + 1)*sizeof(double));

	/*----------------------here are the variables I added----------------------------*/
	/*S is the M x Nt+1 x Nt+1 matrice, storing volatility values*/
	S = (double***)calloc(M, sizeof(double**));
	if (S == NULL)
		return MEMORY_ALLOCATION_FAILURE;
	for (uint j = 0; j < M; j++) /*for each price grid element we generate a markov chain-resided V matrice, size Nt+1 to Nt+1*/
	{
		S[j] = (double**)calloc(M, sizeof(double*));
		if (S[j] == NULL)
			return MEMORY_ALLOCATION_FAILURE;
		for (i = 0; i<Nt + 1; i++)
		{
			S[j][i] = (double *)calloc(Nt + 1, sizeof(double));
			if (S[j][i] == NULL)
				return MEMORY_ALLOCATION_FAILURE;
		}
	}
	return OK;

	ba_log_prices = (double *)calloc(M, sizeof(double));
	if (ba_log_prices == NULL)
	{
		return MEMORY_ALLOCATION_FAILURE;
	}

	fftfreqs = (double *)calloc(M, sizeof(double));
	if (fftfreqs == NULL)
	{
		return MEMORY_ALLOCATION_FAILURE;
	}
}

static void free_memory(uint Nt, uint N, uint M)
{
	uint i;

	for (i = 0; i<Nt + 1; i++)
		free(V[i]);
	free(V);
	for (i = 0; i<Nt + 1; i++)
		free(pu_f[i]);
	free(pu_f);

	for (i = 0; i<Nt + 1; i++)
		free(pd_f[i]);
	free(pd_f);

	for (i = 0; i<Nt + 1; i++)
		free(f_down[i]);
	free(f_down);

	for (i = 0; i<Nt + 1; i++)
		free(f_up[i]);
	free(f_up);

	for (i = 0; i<N + 1; i++)
		free(P_old[i]);
	free(P_old);

	for (i = 0; i<N + 1; i++)
		free(P_new[i]);
	free(P_new);

	for (uint j = 0; j < M; j++) /*for each price grid element we generate a markov chain-resided V matrice, size Nt+1 to Nt+1*/
	{
		for (i = 0; i<Nt + 1; i++)
		{
			free(S[j][i]);
		}
		free(S[j]);
	}
	free(S);

	free(ba_log_prices);
	return;
}

static double compute_f(double r, double omega)
{
	return 2.*sqrt(r) / omega;
}

static double compute_v(double R, double omega)
{
	double val;

	val = SQR(R)*SQR(omega) / 4.;
	if (R>0.)
		val = SQR(R)*SQR(omega) / 4.;
	else
		val = 0.0;
	return val;
}

static double compute_S(double Y, double rv, double omega, double rho)
{
	double val;

	val = exp(Y)*exp(rho*rv / omega);

	return val;
}

static int tree_v(double tt, double v0, double kappa, double theta, double omega, int Nt)
{
	int i, j;
	int z; /*a variable for k_u or k_d, to add to k on n+1 step*/
	double Ru, Rd; /*stores k_u(n,k) and k_d(n,k), respectively*/
	double mu_r, v_curr;
	double dt, sqrt_dt;

	/*Fixed tree for R=f*/
	f[0][0] = compute_f(v0, omega);

	dt = tt / (double)Nt;
	sqrt_dt = sqrt(dt);

	V[0][0] = compute_v(f[0][0], omega);
	f[1][0] = f[0][0] - sqrt_dt;
	f[1][1] = f[0][0] + sqrt_dt;
	V[1][0] = compute_v(f[1][0], omega);
	V[1][1] = compute_v(f[1][1], omega);
	for (i = 1; i<Nt; i++)
		for (j = 0; j <= i; j++)
		{
			f[i + 1][j] = f[i][j] - sqrt_dt;
			f[i + 1][j + 1] = f[i][j] + sqrt_dt;
			V[i + 1][j] = compute_v(f[i + 1][j], omega);
			V[i + 1][j + 1] = compute_v(f[i + 1][j + 1], omega);
		}

	/*Evolve tree for f*/
	for (i = 0; i<Nt; i++)
	{
		for (j = 0; j <= i; j++)
		{
			/*Compute mu_f*/
			v_curr = V[i][j];

			mu_r = kappa*(theta - v_curr);

			z = 0;
			while ((V[i][j] + mu_r*dt<V[i + 1][j - z])
				&& (j - z >= 0)) {

				z = z + 1;
			}
			f_down[i][j] = -z;
			Rd = V[i + 1][j - z];

			z = 0;
			while ((V[i][j] + mu_r*dt>V[i + 1][j + z])
				&& (j + z <= i))
			{
				z = z + 1;
			}

			Ru = V[i + 1][j + z];

			f_up[i][j] = z;
			pu_f[i][j] = (V[i][j] + mu_r*dt - Rd) / (Ru - Rd);

			if ((Ru - 1.e-9>V[i + 1][i + 1]) || (j + f_up[i][j]>i + 1))
			{
				pu_f[i][j] = 1;

				f_up[i][j] = i + 1 - j;
				f_down[i][j] = i - j;
			}

			if ((Rd + 1.e-9<V[i + 1][0]) || (j + f_down[i][j]<0))
			{
				pu_f[i][j] = 0.;
				f_up[i][j] = 1 - j;
				f_down[i][j] = 0 - j;
			}
			pd_f[i][j] = 1. - pu_f[i][j];

		}
	}

	return 1;
}

static int fftfreq(uint n, double d) 
{
	/*
	Return the Discrete Fourier Transform sample frequencies as in fftfreq in numpy
	*/
	double val = 1.0 / (n * d);
	uint N = (n - 1) / 2 + 1;
	int i = 0;
	for (i = 0; i < n-1; i++)
	{
		fftfreqs[i] = val * i;
	}
	for (i = n-1; i < n; i++)
	{
		fftfreqs[i] = val * (-(n / 2) + i);
	}
	return 1;
}

static int compute_price(double tt, double H, double K, double r_premia, double v0, double kappa, double theta, double sigma, double rho, 
	double L, uint M, uint Nt) /*fill matrices with initial conditions and creates basic arrays*/
{
	/*Variables*/
	uint i, j, n, k;
	double r; /*continuous rate*/
	double min_log_price, max_log_price;
	double ds, dt; /*price and time discretization steps*/
	double rho_hat; /*parameter after substitution*/
	double q, factor; /*pde parameters*/

	if (2.0 * kappa * theta < pow(sigma, 2))
		return 1; /*Novikov condition not satisfied*/
	/*Body*/
	r = log(1 + r_premia / 100);
	/*spacial variable. Price space construction*/
	min_log_price = L * log(0.5);
	max_log_price = L * log(2.0);
	ds = (max_log_price - min_log_price) / double(M);
	
	for (j = 1; j < M; j++)
	{
		ba_log_prices[j] = min_log_price + j*ds;
	}
	dt = tt / double(Nt);

	/*building voltree*/
	tree_v(tt, v0, kappa, theta, sigma, Nt);
	
	/*fft frequences we'll need in every vertice of a tree*/
	fftfreq(M, ds);
	rho_hat = sqrt(1.0 - pow(rho, 2.0));
	q = 1.0 / dt + r;
	factor = pow(q*dt, -1.0);

	/*filling F matrice by initial (in time T) conditions*/
	for (j = 0; j < M; j++)
		for (n = 0; n < Nt + k; n++)
			for (k = 0; k < Nt + 1; k++)
			{
				if
				F[j][n][k] = 0;
			}

	return OK;
}
int main()
{
	/*Option parameters*/
	double tt = 1.0;
	double H = 90.0;
	double K = 100.0;
	double r_premia = 10;

	/*Heston model parameters*/
	double V0 = 0.1; /* initial volatility */
	double kappa = 2.0; /*heston parameter, mean reversion*/
	double theta = 0.1; /*heston parameter, long-run variance*/
	double sigma = 0.2; /*heston parameter, volatility of variance*/
	double omega = sigma; /*sigma is used everywhere, omega - in the variance tree*/
	double rho = 0.5; /*heston parameter, correlation*/
	/*method parameters*/
	uint Nt = 100; /*number of time steps*/
	uint M = 512; /*space grid. should be a power of 2*/
	uint L = 3; /*scaling coefficient*/

	memory_allocation(Nt, M, M);
	free_memory(Nt, M, M);
//	printf("Price %.6f  %.6f\n", *price, *delta);
	printf("Hello World!");
	return OK;
}