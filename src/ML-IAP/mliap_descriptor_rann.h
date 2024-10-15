#ifndef LMP_MLIAP_DESCRIPTOR_RANN_H
#define LMP_MLIAP_DESCRIPTOR_RANN_H

#include "mliap_descriptor.h"

namespace LAMMPS_NS{

class MLIAPDescriptorRANN : public MLIAPDescriptor{
public: 
   MLIAPDescriptorRANN(LAMMPS * , char *);
   ~MLIAPDescriptorRANN() override;
   void compute_descriptors (class MLIAPData *) override;
   void compute_forces (class MLIAPData *) override;
   void compute_force_gradients (class MLIAPData *) override;
   void compute_descriptor_gradients (class MLIAPData *) override;
   void init() override;
   double memory_usage() override;  

   void read_paramfile(char *);

  double rcutfac;
  int allocated = 0;
  int max_num = 0;
  char *ctilde_file;

 protected:
  int re;
  int rc;
  int * alpha;
  int dr;
  int o;
  int n;

};

};

#endif