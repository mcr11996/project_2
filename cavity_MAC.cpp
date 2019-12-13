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
		// Solve for F and G
		for (i=1; i<=(n-2); i++){
			for (j=1; j<=(n-1); j++){
				F[i][j] = u[i][j]+dt*((u[i+1][j]-2*u[i][j]+u[i-1][j])/(Re*dx*dx)+
						(u[i][j+1]-2*u[i][j]+u[i][j-1])/(Re*dy*dy)-
						((u[i][j]+u[i+1][j])*(u[i][j]+u[i+1][j])-(u[i][j]+u[i-1][j])*(u[i][j]+u[i-1][j]))/(4*dx)-
						((u[i][j+1]+u[i][j])*(v[i+1][j]+v[i][j])-(u[i][j]+u[i][j-1])*(v[i+1][j-1]+v[i][j-1]))/(4*dy));
			}
		}


		for (i=1; i<=(n-1); i++){
			for (j=1; j<=(n-2); j++){
				G[i][j] = v[i][j]+dt*((v[i+1][j]-2*v[i][j]+v[i-1][j])/(Re*dx*dx)+
						(v[i][j+1]-2*v[i][j]+v[i][j-1])/(Re*dy*dy)-
						((v[i][j]+v[i+1][j])*(u[i][j]+u[i][j+1])-(v[i][j]+v[i-1][j])*(u[i-1][j]+u[i-1][j+1]))/(4*dx)-
						((v[i][j]+v[i][j+1])*(v[i][j]*v[i][j+1])-(v[i][j]+v[i][j-1])*(v[i][j]+v[i][j-1]))/(4*dy));
			}
		}

        // Solve for Pressure
		while (itr < maxitr){
			tempP[i][j] = p[i][j];
			for (i=1; i<=(n-1); i++){
				p[i][0] = p[i][1];
				p[i][n] = p[i][n-1];
			}
			for (j=0; j<=n; j++){
				p[0][j] = p[1][j];
				p[n][j] = p[n-1][j];
			}
			for (i=1; i<=(n-1); i++){
				for (j=1; j<=(n-1); j++){
					p[i][j] = (dy*dy*(p[i-1][j]+p[i+1][j])+
							dx*dx*(p[i][j-1]+p[i][j+1])-
							dx*dx*dy*dy*((F[i][j]-F[i-1][j])/dx+
							(G[i][j]-G[i][j-1])/dy)/dt)/(2*(dx*dx+dy*dy));
					p[i][j] = tempP[i][j] + (p[i][j]-tempP[i][j]);
				}
			}
			if (fabs(p-tempP) <= tol){
				break;
			}
			itr++;
		}


        // Update u:
		for (i=0; i<=(n-1); i++){
			for (j=0; j<=n; j++){
				u[i][j] = F[i][j]-dt*(p[i+1][j]-p[i][j])/dx;
			}
		}

		// Update v:
		for (i=0; i<=n; i++){
			for (j=0; j<= (n-1); j++){
				v[i][j] = G[i][j]-dt*(p[i][j+1]-p[i][j])/dy;
			}
		}
    
	// Values at normal grid points:
	for (i=0; i<=(n-1); i++){
		for (j=0; j<=(n-1); j++){
			uc[i][j] = 0.5*(u[i][j]+u[i][j+1]);
			vc[i][j] = 0.5*(v[i][j]+v[i+1][j]);
			pc[i][j] = 0.25*(p[i][j]+p[i+1][j]+p[i][j+1]+p[i+1][j+1]);
		}
	}
	
	// Output:
	FILE *pfile;
	pfile = fopen("UVP.plt", "w+t");
	if (pfile != NULL){
		fprintf(pfile, "VARIABLES=\"x\",\"Y\",\"U\",\"V\",\"P\"n");
		fprintf(pfile, "ZONE F=POINT\n");
		fprintf(pfile, "I=%d, J=%d\n", n, n);

		for (j=0; j<n; j++){
			for (i=0; i<n; i++){
				double xpos, ypos;
				xpos = i*dx;
				ypos = j*dy;
				fprintf(pfile, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", xpos, ypos, uc[i][j], vc[i][j], pc[i][j]);
			}
		}
	}
	fclose(pfile);


	return 0;
}

	 


}
