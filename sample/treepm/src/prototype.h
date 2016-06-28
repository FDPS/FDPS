#ifndef __TREEPM_PROTOTYPE_H__
#define __TREEPM_PROTOTYPE_H__

extern "C"{
  void make_directory(char*);
  float wallclock_timing(struct timeval, struct timeval);
  float timing(struct tms, struct tms);
}

#endif
