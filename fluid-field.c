#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// --- Simulation Parameters ---
#define NX 101   // CHANGED: Grid points in X (was 41)
#define NY 31    // CHANGED: Grid points in Y (was 41)
#define NIT 100  // CHANGED: Pressure iterations (was 50)
#define NT 4000  // CHANGED: Total time steps (was 100)
#define T 4.0    // CHANGED: Actual time in seconds (was 0.1)
#define SKIP 20  // CHANGED: How often to write data (was 1)

// --- Physics Parameters ---
double lx = 6.0; // channel length 
double ly = 2.0; // channel height 
double rho = 1.0;
double nu = 0.3; 
double u_0 = 1.0; // inlet velocity 

double dt = 1.0*T/NT;

// Obstacle parameters for vertical plate
double obs_cx = 1.26;   // center x
double obs_cy = 1.0;    // center y  
double obs_w = 0.15;    // width (thickness)
double obs_h = 0.7;     // height (length)

// Obstacle mask array (1 = obstacle, 0 = fluid)
int obstacle[NY][NX];

// Arrays
double u[NY][NX], v[NY][NX], p[NY][NX];
double un[NY][NX], vn[NY][NX], pn[NY][NX];
double b[NY][NX];

// Function to check if point are inside rectangle
int inside_obstacle(double x, double y) {
    double dx = x - obs_cx;
    double dy = y - obs_cy;
    return (fabs(dx) <= obs_w/2.0 && fabs(dy) <= obs_h/2.0);
}

int main() {
    double dx = lx/(NX-1);
    double dy = ly/(NY-1);

    // Build obstacle mask by checking each grid point
    for(int j=0; j<NY; j++) {
        for(int i=0; i<NX; i++) {
            double x = i * dx;
            double y = j * dy;
            obstacle[j][i] = inside_obstacle(x, y);
        }
    }

    // Initialize arrays - start from rest
    for(int j=0; j<NY; j++) {
        for(int i=0; i<NX; i++) {
            u[j][i] = 0.0;
            v[j][i] = 0.0; 
            p[j][i] = 0.0;
            un[j][i] = 0.0; 
            vn[j][i] = 0.0; 
            pn[j][i] = 0.0;
            b[j][i] = 0.0;
        }
    }

    fprintf(stderr, "Starting Channel Flow Simulation...\n");

    // --- MAIN TIME LOOP ---
    for (int n=0; n <= NT; n++) {

        // 1. Compute Source Term (b) for Pressure Equation
        for (int j=1; j < NY-1; j++) {
            for (int i=1; i < NX-1; i++) {
                if (obstacle[j][i]) continue;  // ADDED: Skip obstacle cells
                b[j][i] = rho*(1.0/dt)*((u[j][i+1]-u[j][i-1])/(2.0*dx)+
                                       (v[j+1][i]-v[j-1][i])/(2.0*dy));
            }
        }

        // 2. Pressure Poisson Solver
        for (int it = 0; it < NIT; it++) {
            for(int j=0; j<NY; j++) 
                for(int i=0; i<NX; i++) 
                    pn[j][i] = p[j][i];

            for (int j=1; j < NY-1; j++) {
                for (int i=1; i < NX-1; i++) {
                    if (obstacle[j][i]) continue;  // ADDED: Skip obstacle cells
                    p[j][i] = ((pn[j][i+1]+pn[j][i-1])*dy*dy + 
                               (pn[j+1][i]+pn[j-1][i])*dx*dx -
                               b[j][i]*dx*dx*dy*dy) / (2.0*(dx*dx+dy*dy));
                }
            }

            // Pressure BCs - pressure driven flow
            for (int j=0; j < NY; j++) {
                p[j][0] = 1.0;           // CHANGED: inlet high pressure (drives flow)
                p[j][NX-1] = 0;          // outlet: p = 0
            }
            for (int i=0; i < NX; i++) {
                p[0][i] = p[1][i];       // bottom wall
                p[NY-1][i] = p[NY-2][i]; // top wall
            }
        }

        // Copy u, v to un, vn
        for(int j=0; j<NY; j++) {
            for(int i=0; i<NX; i++) {
                un[j][i] = u[j][i];
                vn[j][i] = v[j][i];
            }
        }

        // 3. Update Velocities
        for (int j=1; j < NY-1; j++) {
            for (int i=1; i < NX-1; i++) {
                if (obstacle[j][i]) continue;  // ADDED: Skip obstacle cells
                    
                u[j][i] = un[j][i] 
                    - un[j][i]*dt/dx*(un[j][i]-un[j][i-1]) 
                    - vn[j][i]*dt/dy*(un[j][i]-un[j-1][i]) 
                    - dt/(2.0*rho*dx)*(p[j][i+1]-p[j][i-1]) 
                    + nu*dt/(dx*dx)*(un[j][i+1]-2.0*un[j][i]+un[j][i-1]) 
                    + nu*dt/(dy*dy)*(un[j+1][i]-2.0*un[j][i]+un[j-1][i]);

                v[j][i] = vn[j][i] 
                    - un[j][i]*dt/dx*(vn[j][i]-vn[j][i-1]) 
                    - vn[j][i]*dt/dy*(vn[j][i]-vn[j-1][i]) 
                    - dt/(2.0*rho*dy)*(p[j+1][i]-p[j-1][i]) 
                    + nu*dt/(dx*dx)*(vn[j][i+1]-2.0*vn[j][i]+vn[j][i-1]) 
                    + nu*dt/(dy*dy)*(vn[j+1][i]-2.0*vn[j][i]+vn[j-1][i]);
            }
        }

        // 4. Velocity Boundary Conditions
        
        //Top and bottom walls: free-slip (not no-slip)
        for (int i=0; i < NX; i++) {
            u[0][i] = u[1][i];      v[0][i] = 0;      // bottom: du/dy=0, v=0
            u[NY-1][i] = u[NY-2][i];  v[NY-1][i] = 0; // top: du/dy=0, v=0
        }
        
        //Inlet - velocity develops from pressure gradient
        for (int j=1; j < NY-1; j++) {
            u[j][0] = u[j][1];  // copy from inside
            v[j][0] = 0;
        }
        
        // Outlet: copy from inside
        for (int j=0; j < NY; j++) {
            u[j][NX-1] = u[j][NX-2];
            v[j][NX-1] = v[j][NX-2];
        }
        
        // Obstacle no-slip boundary condition
        for (int j=0; j < NY; j++) {
            for (int i=0; i < NX; i++) {
                if (obstacle[j][i]) {
                    u[j][i] = 0;
                    v[j][i] = 0;
                }
            }
        }

        // Output data
        if (n % SKIP == 0) {
            for (int j = 0; j < NY; j++) {
                for (int i = 0; i < NX; i++) {
                    printf("%f %f %f %f\n", i*dx, j*dy, u[j][i], v[j][i]);
                }
            }
            printf("\n\n");
        }
    }

    fprintf(stderr, "Simulation Complete.\n");
    return 0;
}
