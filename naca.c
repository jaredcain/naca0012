/* I believe these are the only header files i need to include for the simple flow to find lift and drag */

#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

#define chord       (1.00)
#define uref        (1.00)
#define thickness   (0.12)

double Re = 6.e6;
double aoa = 0. * M_PI / 180.0;
int level = 12;

/* I've seen this in many examples, unsure what it does beyond mu being viscosity */

face vector muv[];

/* This is how i build the airfoil shape, unsure if this is sufficient, as it was mentioned that the trailing edge should be treated differently */

double naca0012(double x, double y) {
  if (x >= 0. && x <= 1.) {
    return sq(y) -sq(5. * thickness * (0.2969 * sqrt(x) - 0.1260 * x - 0.3516 * sq(x) + 0.2843 * cube(x) - 0.1036 * pow(x, 4.)));
   } 
  else {
    return 1.0;
  }
}

/* This is how aoa is applied, works well from what i can tell */

double naca(double x, double y)
{
  double xr = x * cos(aoa) - y * sin(aoa);
  double yr = x * sin(aoa) + y * cos(aoa);
  return naca0012(xr, yr);
}

/* main function where i set up being able to choose AoA, have the logfile, set the domain size and number of cells, unsure again what mu = muv does, this is the theory side that I'm lacking on */

int main(int argc, char *argv[])
{
  if (argc > 1) {
    aoa = atof(argv[1]) * M_PI / 180.0;
  }
  
  char log_filename[50];
  snprintf(log_filename, sizeof(log_filename), "%.0f.log", aoa * 180.0 / M_PI);
  freopen(log_filename, "w", stderr);
  
  L0 = 16.;
  origin (-4, -L0/2.);
  N = 1 << level;
  mu = muv;
  run(); 
}

/* I understand this is viscosity, unsure if it is the only think needed, and if it is done correctly */

event properties (i++)
{
  foreach_face()
    muv.x[] = (uref)*(chord)/(Re)*fm.x[];
}

/* BC's, I think some may  be redundant and/or unnecessary */

u.n[left]  = dirichlet(1);
p[left]    = neumann(0);
pf[left]   = neumann(0);

u.n[right] = neumann(0);
p[right]   = dirichlet(0);
pf[right]  = dirichlet(0);

u.n[embed] = dirichlet(0);
u.t[embed] = dirichlet(0);

/* Initialization, this embeds the airfoil shape? think u.x=1 may be redundant based on BCs */

event init (t = 0) {
  solid(cs, fs, naca(x, y));
  foreach()
    u.x[] = 1.;
}

/* Fairly certain all is correct here, but I'm still having some issues with negative drag coefficients */

event logfile (i += 10) {
  scalar omega[];
  vorticity (u, omega);
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x + Fmu.x) / (0.5 * sq((uref)) * (chord));
  double CL = (Fp.y + Fmu.y) / (0.5 * sq((uref)) * (chord));
  fprintf(stderr, "%.4d\t%.4e\t%.4e\t%.4e\n", i, t, CD, CL);
}

event movies (t += 0.05, t <= 30)
{
  scalar p[], m[];
  view(quat = {0.000, 0.000, 0.000, 1.000}, fov = 30, near = 0.01, far = 1000, tx = -0.10, ty = 0.00, tz = -0.35, width = 1024, height = 512);
  box();
  squares(color = "p", linear=true); 
  draw_vof(c = "cs", lw = 1, lc = {0, 0, 0});
  char video[50];
  snprintf(video, sizeof(video), "%.0f.mp4", aoa * 180.0 / M_PI);
  save(video);
}

/* Adaptive Mesh, not sure where minlevel should be, if tolerances are low enough and/or too low */

event adapt (i++) {
  adapt_wavelet ({cs,u.x,u.y}, (double[]){1e-2,5e-3,5e-3}, maxlevel=level, minlevel=4);
}
