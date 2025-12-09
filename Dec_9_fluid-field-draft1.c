#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// --- Simulation Parameters ---
#define NX 328     // Grid points in X
#define NY 82    // Grid points in Y
#define NIT 50   // Pressure iterations per time step
#define NT 120   // Total time steps
#define T 0.1    // Actual time in seconds
#define SKIP 1   // How often to write data
#define VI 0.2 // The starting speed of the pipe
#define XI 25 // Starting X coord for the box
#define XF 40 // Ending X coord for the box

// --- Physics Parameters ---
// Reynolds Number Re = (U*L)/nu. 
// High Re = Turbulent, Low Re = Viscous
// Time   Re  Viscosity (nu) Fluid Analogy
// ----   --  -------------- -------------
// 0.06s  0.3   6.67 m^2/s    Molten Glass
// 0.15s  1     2.0  m^2/s    Glycerin (cold)
//  1.5s  10    0.2  m^2/s    SAE 40 Oil
//  7.5s  50    0.04 m^2/s    Vegetable Oil
//
// Keep Re low for this simple solver code stability.

double lx = 16.0; // physical length in X direction
double ly = 4.0; // physical length in Y direction
double rho = 1.0;
double nu = 0.3; // extremely viscous 
double dt = 1.0*T/NT;

// Arrays
double u[NY][NX], v[NY][NX], p[NY][NX];  // Main variables
double un[NY][NX], vn[NY][NX], pn[NY][NX]; // Previous time step buffers
double b[NY][NX];              // Source term for pressure equation

int main() {
  double dx = lx/(NX-1);
  double dy = ly/(NY-1);

  // Initialize arrays to zero
  for(int j=0; j<NY; j++) {
    for(int i=0; i<NX; i++) {
      u[j][i] = VI; v[j][i] = 0.0; p[j][i] = 0.0;
      un[j][i] = 0; vn[j][i] = 0.0; pn[j][i] = 0.0;
      b[j][i] = 0.0;
    }
  }

  printf("Starting Simulation...\n");

  // --- MAIN TIME LOOP ---
  for (int n=0; n <= NT; n++) {
    
    // 1. Compute Source Term (b) for Pressure Equation
    // This represents the divergence of the tentative velocity
    for (int j=1; j < NY-1; j++) {
      for (int i=1; i < NX-1; i++) {
        b[j][i] = (rho*(1.0/dt)*((u[j][i+1]-u[j][i-1])/(2.0*dx)+
                               (v[j+1][i]-v[j-1][i])/(2.0*dy)));
        
     // Could add non-linear terms to source here but whatev's

      }
    }

    // 2. Pressure Poisson Solver
    // We iterate NIT times to smooth out the pressure field
    for (int it = 0; it < NIT; it++) {
      // Copy p to pn
      for(int j=0; j<NY; j++) 
        for(int i=0; i<NX; i++) 
          pn[j][i] = p[j][i];

      for (int j=1; j < NY-1; j++) {
        for (int i=1; i < NX-1; i++) {
          p[j][i] = ((pn[j][i+1]+pn[j][i-1])*dy*dy + 
                 (pn[j+1][i]+pn[j-1][i])*dx*dx -
                 b[j][i]*dx*dx*dy*dy)/(2.0*(dx*dx+dy*dy));
        }
      }

      // Pressure Boundary Conditions (Neumann: dp/dn = 0)
      // Copy neighbor value to wall
      for (int j=0; j < NY; j++) {
        p[j][0] = p[j][1];      // Left
        p[j][NX-1] = p[j][NX-2];  // Right
      }
      for (int i=0; i < NX; i++) {
        p[0][i] = p[1][i];      // Bottom
        p[NY-1][i] = p[NY-2][i];  // Top
      }
    }

    // Copy current u, v to un, vn
    for(int j=0; j<NY; j++) {
      for(int i=0; i<NX; i++) {
        un[j][i] = u[j][i];
        vn[j][i] = v[j][i];
      }
    }

    // 3. Update Velocities (Advection + Diffusion + Pressure Correction)
    for (int j=1; j < NY-1; j++) {
      for (int i=1; i < NX-1; i++) {
        
        // Navier-Stokes Momentum Equations
        
        // u-momentum
        u[j][i] = un[j][i] 
            - un[j][i]*dt/dx*(un[j][i]-un[j][i-1]) 
            - vn[j][i]*dt/dy*(un[j][i]-un[j-1][i]) 
            - dt/(2.0*rho*dx)*(p[j][i+1]-p[j][i-1]) 
            + nu*dt/(dx*dx)*(un[j][i+1]-2.0*un[j][i]+un[j][i-1]) 
            + nu*dt/(dy*dy)*(un[j+1][i]-2.0*un[j][i]+un[j-1][i]);

        // v-momentum
        v[j][i] = vn[j][i] 
            - un[j][i]*dt/dx*(vn[j][i]-vn[j][i-1]) 
            - vn[j][i]*dt/dy*(vn[j][i]-vn[j-1][i]) 
            - dt/(2.0*rho*dy)*(p[j+1][i]-p[j-1][i]) 
            + nu*dt/(dx*dx)*(vn[j][i+1]-2.0*vn[j][i]+vn[j][i-1]) 
            + nu*dt/(dy*dy)*(vn[j+1][i]-2.0*vn[j][i]+vn[j-1][i]);
      }
    }

    // 4. Velocity Boundary Conditions
    // Stationary Walls
    
			for (int j=0; j < NY; j++) {
      u[j][0] = VI; u[j][NX] = 0;//left wall ; right wall (x velocity) 
      v[j][0] = 0; v[j][NX] = 0;//left wall : right wall (y velocity)
    }
    for (int i=0; i < NX; i++) {
      u[0][i] = 0; v[0][i] = 0; //bottom x velocity; bottom y velocity 
      u[NX][i] = 0; v[NX][i] = 0; //bottom x velocity; bottom y velocity 

      u[NY-1][i] = 0;//top wall x velocity  
      v[NY-1][i] = 0; //top wall y velocity 

      u[0][i] = 0;//top wall x velocity  
      v[0][i] = 0; //top wall y velocity 
    }
    //BLOCKER
		
		// Set the inside to 0
		for (int i = NY/4+1; i< 3*NY/4-1; i++){
			for (int k = XI+1; k < XF-1;k++){
				u[i][k] = 0;
				v[i][k] =0;
			}
		}

		for (int i=NY/4; i < (3*NY/4); i++){ //front vertical portion
			u[i][XI] = 0;
			v[i][XI] = v[i][XI-1];
		}	
	
		for (int i=NY/4; i < (3*NY/4); i++){ //rear vertical portion
			u[i][XF] = 0;
			v[i][XF] = v[i][XI+1];
		}	
			
		for (int i=XI; i < XF; i++){ //Top portion
			u[3*NY/4][i] = u[3*NY/4][i-1];
			v[3*NY/4][i] = 0;
		}	
		

		for (int i=XI; i < XF; i++){// Bottom Portion
			u[NY/4][i] = u[NY/4-1][i];
			v[NY/4][i] = 0;
		}	
    if (n % SKIP == 0) {
				for (int j = 0; j < NY; j++) {
					  for (int i = 0; i < NX; i++) {
				   printf("%f %f %f %f\n",i*dx,j*dy,u[j][i],v[j][i]);
				  }
				}
				printf("\n\n");
    }
  }

  fprintf(stderr,"Simulation Complete.\n");
  return 0;
}
