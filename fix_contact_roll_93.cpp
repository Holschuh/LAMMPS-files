#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_contact_roll_93.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

const char syntax[] =
  "fix contact/roll93: syntax is "
  "fix ID group-ID contact/roll93 "
  "<cx|v_> <cy|v_> <cz|v_> <ux> <uy> <uz> "
  "<radius> <epsilon> <sigma> <cutoff>";

/* ---------------------------------------------------------------------- */

FixContactRoll93::FixContactRoll93(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  char *endptr;

  if (narg < 13)
    error->all(FLERR,syntax);

  cxvar = cyvar = czvar = NULL;

  // --- cx ---
  if (strstr(arg[3], "v_") == arg[3]) {
    cxvar = strdup(arg[3]);
  } else {
    cx = strtod(arg[3], &endptr);
    if (endptr == arg[3]) error->all(FLERR,syntax);
  }

  // --- cy ---
  if (strstr(arg[4], "v_") == arg[4]) {
    cyvar = strdup(arg[4]);
  } else {
    cy = strtod(arg[4], &endptr);
    if (endptr == arg[4]) error->all(FLERR,syntax);
  }

  // --- cz ---
  if (strstr(arg[5], "v_") == arg[5]) {
    czvar = strdup(arg[5]);
  } else {
    cz = strtod(arg[5], &endptr);
    if (endptr == arg[5]) error->all(FLERR,syntax);
  }

  // Richtungsvektor (nur numerisch)
  ux = strtod(arg[6], &endptr);
  uy = strtod(arg[7], &endptr);
  uz = strtod(arg[8], &endptr);

  if (endptr == arg[6] || endptr == arg[7] || endptr == arg[8])
    error->all(FLERR,syntax);

  // Normieren
  unorm = sqrt(ux*ux + uy*uy + uz*uz);
  if (unorm == 0.0)
    error->all(FLERR,"fix contact/roll93: zero direction vector");

  ux /= unorm;
  uy /= unorm;
  uz /= unorm;

  // Parameter
  radius  = strtod(arg[9], &endptr);
  epsilon = strtod(arg[10], &endptr);
  sigma   = strtod(arg[11], &endptr);
  cutoff  = strtod(arg[12], &endptr);

  if (endptr == arg[9] || endptr == arg[10] ||
      endptr == arg[11] || endptr == arg[12])
    error->all(FLERR,syntax);

  // Precompute
  radius_sq = radius*radius;
  radius_plus_cutoff_sq = (radius+cutoff)*(radius+cutoff);

  coeff1 = 18.0/15.0 * epsilon * pow(sigma,9.0);
  coeff2 = 3.0 * epsilon * pow(sigma,3.0);
  coeff3 = 2.0/15.0 * epsilon * pow(sigma,9.0);
  coeff4 = 1.0 * epsilon * pow(sigma,3.0);

  double r3inv = 1.0/(cutoff*cutoff*cutoff);
  double r9inv = r3inv*r3inv*r3inv;
  offset = coeff3*r9inv - coeff4*r3inv;

  energy_global_flag = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 4;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  force_flag = 0;
  for (int i = 0; i < 5; i++) {
    froll[i] = 0.0;
    froll_loc[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

FixContactRoll93::~FixContactRoll93()
{
  if (cxvar) free(cxvar);
  if (cyvar) free(cyvar);
  if (czvar) free(czvar);
}

/* ---------------------------------------------------------------------- */

int FixContactRoll93::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixContactRoll93::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // Variablen auswerten
  if (cxvar) cx = input->variable->compute_equal(cxvar);
  if (cyvar) cy = input->variable->compute_equal(cyvar);
  if (czvar) cz = input->variable->compute_equal(czvar);

  for (int i = 0; i < 5; i++) froll_loc[i] = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      double rx = x[i][0] - cx;
      double ry = x[i][1] - cy;
      double rz = x[i][2] - cz;

      double delta[3] = {rx, ry, rz};
      domain->minimum_image(FLERR, delta);
      rx = delta[0];
      ry = delta[1];
      rz = delta[2];

      // Projektion
      double dot = rx*ux + ry*uy + rz*uz;

      double dx = rx - dot*ux;
      double dy = ry - dot*uy;
      double dz = rz - dot*uz;

      double rsq = dx*dx + dy*dy + dz*dz;

      if (rsq < radius_sq)
        error->one(FLERR,"fix contact/roll93: atom inside cylinder");

      if (rsq < radius_plus_cutoff_sq) {

        double r = sqrt(rsq);
        double rinv = 1.0/(r - radius);
        double r3inv = rinv*rinv*rinv;
        double r9inv = r3inv*r3inv*r3inv;

        double df = (coeff1*r9inv - coeff2*r3inv) * rinv;
        double e  = (coeff3*r9inv - coeff4*r3inv) - offset;

        double fx = df * dx / r;
        double fy = df * dy / r;
        double fz = df * dz / r;

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        froll_loc[0] += e;
        froll_loc[1] -= fx;
        froll_loc[2] -= fy;
        froll_loc[3] -= fz;
        froll_loc[4] += 1.0;
      }
    }
  }

  force_flag = 0;
}




/* ---------------------------------------------------------------------- */

void FixContactRoll93::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    error->all(FLERR,"fix contact/roll93: RESPA not supported.");
}

/* ---------------------------------------------------------------------- */

void FixContactRoll93::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixContactRoll93::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixContactRoll93::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

double FixContactRoll93::compute_scalar()
{
  if (force_flag == 0) {
    MPI_Allreduce(froll_loc, froll, 5, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return froll[0];
}

/* ---------------------------------------------------------------------- */

double FixContactRoll93::compute_vector(int n)
{
  if (force_flag == 0) {
    MPI_Allreduce(froll_loc, froll, 5, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return froll[n+1];
}






