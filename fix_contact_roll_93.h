#ifdef FIX_CLASS

FixStyle(contact/roll,FixContactRoll93)

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
  // Variablen für dynamische Zentren (v_*)
  char *cxvar, *cyvar, *czvar;

  // Zylinderzentrum
  double cx, cy, cz;

  // Richtungsvektor des Zylinders
  double ux, uy, uz;
  double unorm;

  // Potentialparameter
  double radius, epsilon, sigma, cutoff;

  // vorab berechnete Werte für Effizienz
  double radius_sq, radius_plus_cutoff_sq;
  double coeff1, coeff2, coeff3, coeff4, offset;

  // interne Flags und Kräfte
  int force_flag;
  double froll[5], froll_loc[5];
};

}

#endif
#endif
