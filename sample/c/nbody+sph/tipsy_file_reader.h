#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int read_ptcl_number(char file_name[]);
void read_ptcl_data(char file_name[], int n_ptcl, float *mass, float *pos, float *vel);

#ifdef __cplusplus
} // END of extern "C"
#endif
