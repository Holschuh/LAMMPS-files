#ifdef FIX_CLASS

FixStyle(contact/roll93,FixContactRoll93)

#else

#ifndef LMP_FIX_CONTACT_ROLL_93_H
#define LMP_FIX_CONTACT_ROLL_93_H

#include "fix.h"

namespace LAMMPS_NS {

class FixContactRoll93 : public Fix {
 public:
  FixContactRoll93(class LAMMPS *, int, char **);
  ~FixContactRoll93();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void min_post_force(int);

  double compute_scalar();
  double compute_vector(int);

 private:
  // Variablen (v_*)
  char *cxvar, *cyvar, *czvar;

  // Zentrum
  double cx, cy, cz;

  // Zylinderachse
  double ux, uy, uz;
  double unorm;

  // LJ Parameter
  double radius, epsilon, sigma, cutoff;

  double radius_sq, radius_plus_cutoff_sq;
  double coeff1, coeff2, coeff3, coeff4, offset;

  int force_flag;
  double froll[5], froll_loc[5];
};

}

#endif
#endif
