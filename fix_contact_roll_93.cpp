#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_contact_roll.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

const char syntax[] =
  "fix contact/roll: Illegal fix command, syntax is 'fix ID group-ID "
  "contact/roll <cx> <cy> <cz> <ux> <uy> <uz> <radius> <epsilon> <sigma> <cutoff>'";

/* ---------------------------------------------------------------------- */

FixContactRoll::FixContactRoll(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  int iarg;
  char *endptr;

  if (narg != 13)
    error->all(FLERR,syntax);

  cxvar = cyvar = czvar = NULL;

  // Zentrum
  cx = strtod(arg[3], &endptr);
  cy = strtod(arg[4], &endptr);
  cz = strtod(arg[5], &endptr);

  if (endptr == arg[3] || endptr == arg[4] || endptr == arg[5])
    error->all(FLERR,syntax);

  // Richtungsvektor
  ux = strtod(arg[6], &endptr);
  uy = strtod(arg[7], &endptr);
  uz = strtod(arg[8], &endptr);

  if (endptr == arg[6] || endptr == arg[7] || endptr == arg[8])
    error->all(FLERR,syntax);

  // Normieren
  unorm = sqrt(ux*ux + uy*uy + uz*uz);
  if (unorm == 0.0)
    error->all(FLERR,"fix contact/roll: zero direction vector");

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

  coeff1 = 48.0 * epsilon * pow(sigma,12.0);
  coeff2 = 24.0 * epsilon * pow(sigma,6.0);
  coeff3 = 4.0 * epsilon * pow(sigma,12.0);
  coeff4 = 4.0 * epsilon * pow(sigma,6.0);

  double r2inv = 1.0/(cutoff*cutoff);
  double r6inv = r2inv*r2inv*r2inv;
  offset = r6inv*(coeff3*r6inv - coeff4);

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

FixContactRoll::~FixContactRoll() {}

/* ---------------------------------------------------------------------- */

int FixContactRoll::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixContactRoll::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    error->all(FLERR,"fix contact/roll: RESPA not supported.");
}

/* ---------------------------------------------------------------------- */

void FixContactRoll::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixContactRoll::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixContactRoll::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < 5; i++)
    froll_loc[i] = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      double rx = x[i][0] - cx;
      double ry = x[i][1] - cy;
      double rz = x[i][2] - cz;

      domain->minimum_image(rx, ry, rz);

      // Projektion auf Achse
      double dot = rx*ux + ry*uy + rz*uz;

      double rpx = dot * ux;
      double rpy = dot * uy;
      double rpz = dot * uz;

      // radialer Abstand zum Zylinder
      double dx = rx - rpx;
      double dy = ry - rpy;
      double dz = rz - rpz;

      double rsq = dx*dx + dy*dy + dz*dz;

      if (rsq < radius_sq) {
        error->one(FLERR,"fix contact/roll: atom inside cylinder");
      }

      if (rsq < radius_plus_cutoff_sq) {

        double r = sqrt(rsq);
        double rinv = 1.0/(r - radius);
        double r2inv = rinv*rinv;
        double r6inv = r2inv*r2inv*r2inv;

        double df = r6inv*(coeff1*r6inv - coeff2) * rinv;
        double e  = r6inv*(coeff3*r6inv - coeff4) - offset;

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

void FixContactRoll::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

double FixContactRoll::compute_scalar()
{
  if (force_flag == 0) {
    MPI_Allreduce(froll_loc, froll, 5, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return froll[0];
}

/* ---------------------------------------------------------------------- */

double FixContactRoll::compute_vector(int n)
{
  if (force_flag == 0) {
    MPI_Allreduce(froll_loc, froll, 5, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return froll[n+1];
}
