#include <iostream>
#include <cmath>

int main() {
	int n = 100;
	double u[n][n+1], F[n][n+1], uc[n][n];
	double v[n+1][n], G[n+1][n], vc[n][n];
	double p[n+1][n+1], tempP[n+1][n+1], pc[n][n];
	int i, j, itr=1, maxitr=500000;
	double dx=1.0/(n-1), dy=1.0/(n-1), dt=0.01, error=1.0, Re=100.0;
	double tol = 0.00001;

	// Initialization
	for (i=0; i<=(n-1); i++){
		for (j=0; j<=n; j++){
			u[i][j] = 0.0;
			u[i][n] = 1.0;
			u[i][n-1] = 1.0;
		}
	}

	for (i=0; i<=(n); i++){
		for (j=0; j<=(n-1); j++){
			v[i][j] = 0.0;
		}
	}

	for (i=0; i<=n; i++){
		for (j=0; j<=n; j++){
			p[i][j] = 1.0;
		}
	}


	while (error > 0.00000001){
	    // Boundary Conditions
		for (j=1; j<=(n-1); j++){
			F[0][j] = 0.0;
			F[n-1][j] = 0.0;
			u[0][j] = 0.0;
			u[n-1][j] = 0.0;
		}

		for (i=0; i<=(n-1); i++){
			F[i][0] = -F[i][1];
			F[i][n] = 2 - F[i][n-1];
			u[i][0] = -u[i][1];
			u[i][n] = 2 - u[i][n-1];
		}

		for (j=1; j<=(n-2); j++){
			G[0][j] = - G[1][j];
			G[n][j] = -G[n-2][j];
			v[0][j] = - v[1][j];
			v[n][j] = -v[n-2][j];
		}

		for (i=0; i<=n; i++){
			G[i][0] = 0.0;
			G[i][n-1] = 0.0;
			v[i][0] = 0.0;
			v[i][n-1] = 0.0;

		}

	return 0;
}

}
