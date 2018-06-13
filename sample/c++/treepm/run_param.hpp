#ifndef __RUN_PARAM_HPP__
#define __RUN_PARAM_HPP__

#include "constants.hpp"

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

#include <cstdio>
#include <cstdlib>
#include <sys/time.h>
#include <sys/times.h>
#include <particle_simulator.hpp>
#include "cosmology.hpp"

#include "prototype.h"

class run_param {
public:
    enum {
        SANTABARBARA = 0,
        READ_FILE = 1,
        RANDOM = 2,
    };
    PS::S32 mode;
    PS::F64 theta;
    PS::S32 step;
    PS::S32 mpi_nproc, mpi_rank;
  
    PS::S64 npart_total;
    PS::S64 npart_local;

    COSM::REAL zend;
    COSM::REAL tnow;
    COSM::REAL anow;
    COSM::REAL znow;
    COSM::REAL hnow;

    COSM::REAL tend;

    PS::F64 lunit, munit, tunit;

    COSM::cosmology cosm;

    char model_name[64];
    FILE *diag_fp;
    PS::S32 noutput, output_indx;
    COSM::REAL *output_timing;

    /*
    void input_modelname(char *filename) {
        param_fp = fopen(filename, "r");
        if (param_fp == NULL) {
            fprintf(stderr, "File %s not found in input_params.\n", filename);
            exit(EXIT_FAILURE);
        }

        int ret = fscanf(param_fp,"%s",model_name);
        if(ret == EOF) 
          fprintf(stderr,"Input error of a model name in input_modelname()\n");
        
        fclose(param_fp);
    }
    */

    void input_params(FILE * param_fp) {
        static char diag_filename[256];
        int ret;

        ret = fscanf(param_fp, "%s", model_name);
        if(ret == EOF) fprintf(stderr,"Input error of a model name in input_params()\n");

        sprintf(diag_filename,"%s.diag", model_name);

        /* Open the diag_file */
        if (mpi_rank == 0) open_diag_file(diag_filename);

        /* # of output timing */
        ret = fscanf(param_fp,"%d",&noutput);
        if (ret == EOF) fprintf(stderr,"Input error of noutout in input_params()\n");

        output_timing = (COSM::REAL *)malloc(sizeof(COSM::REAL)*noutput);
        output_indx = 0;

        COSM::REAL prev_output_timing=0.0;

        for (PS::S32 iout=0;iout<noutput;iout++) {
            int ret = fscanf(param_fp, "%f", &output_timing[iout]);
            if (ret == EOF) 
                fprintf(stderr,"Input error of output_timing in input_params()\n");

            if (iout>0 && output_timing[iout] > prev_output_timing) {
                fprintf(stderr,"Output timing must be in reducing order.\n");
                MPI_Finalize();
                std::exit(EXIT_FAILURE);
            }
            if (output_timing[iout] > znow) {
                output_indx = iout+1;
            }
            prev_output_timing = output_timing[iout];
        }

        tend = cosm.ztotime(output_timing[noutput-1]);

        if (output_indx == noutput || tnow >= tend) {
            if (mpi_rank == 0) {
                fprintf(stderr,"This run already finished.\n");
                fprintf(stderr,"No outputs requested.\n");
            }
        }
    }
    
#if 0    
    void input_params(char *filename) {
        static char diag_filename[256];
        int ret;

        param_fp = fopen(filename, "r");
        if (param_fp == NULL) {
            fprintf(stderr, "File %s not found in input_params.\n", filename);
            exit(EXIT_FAILURE);
        }

        ret = fscanf(param_fp, "%s", model_name);
        if (ret == EOF) fprintf(stderr, "Input error of a model name in input_params()\n");

        sprintf(diag_filename,"%s.diag", model_name);

        /* Open the diag_file */
        if (mpi_rank == 0) open_diag_file(diag_filename);

        /* # of output timing */
        ret = fscanf(param_fp,"%d",&noutput);
        if (ret == EOF) fprintf(stderr, "Input error of noutout in input_params()\n");

        output_timing = (COSM::REAL *)malloc(sizeof(COSM::REAL)*noutput);
        output_indx = 0;

        COSM::REAL prev_output_timing=0.0;

        for (PS::S32 iout=0; iout<noutput; iout++) {
            int ret = fscanf(param_fp, "%f", &output_timing[iout]);
            if (ret == EOF) fprintf(stderr,"Input error of output_timing in input_params()\n");

            if (iout>0 && output_timing[iout] > prev_output_timing) {
                fprintf(stderr, "Output timing must be in reducing order.\n");
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
            if (output_timing[iout] > znow) {
                output_indx = iout+1;
            }
            prev_output_timing = output_timing[iout];
        }

        tend = cosm.ztotime(output_timing[noutput-1]);

        if (output_indx == noutput || tnow >= tend) {
            if (mpi_rank == 0) {
                fprintf(stderr,"This run already finished.\n");
                fprintf(stderr,"No outputs requested.\n");
            }
        }

    }
#endif
    
    void write_header(FILE *fp) {

        if (fp==NULL) {exit(EXIT_FAILURE);}
        int cnt=0;

        cnt += fwrite(&mpi_nproc, sizeof(mpi_nproc),1,fp);
        cnt += fwrite(&mpi_rank,sizeof(mpi_rank),1,fp);
        cnt += fwrite(&npart_local,sizeof(npart_local),1,fp);
        cnt += fwrite(&npart_total,sizeof(npart_total),1,fp);

        cnt += fwrite(&cosm.omegam,sizeof(cosm.omegam),1,fp);
        cnt += fwrite(&cosm.omegab,sizeof(cosm.omegab),1,fp);
        cnt += fwrite(&cosm.omegav,sizeof(cosm.omegav),1,fp);
        cnt += fwrite(&cosm.omeganu,sizeof(cosm.omeganu),1,fp);
        cnt += fwrite(&cosm.hubble,sizeof(cosm.hubble),1,fp);
        cnt += fwrite(&anow, sizeof(anow),1,fp);
        cnt += fwrite(&tnow, sizeof(tnow),1,fp);
        
        cnt += fwrite(&lunit,sizeof(lunit),1,fp);
        cnt += fwrite(&munit,sizeof(munit),1,fp);
        cnt += fwrite(&tunit,sizeof(tunit),1,fp);

    }

    void read_header(FILE *fp) {

        if (fp==NULL) {exit(EXIT_FAILURE);}

        int cnt = 0;
          
        cnt += fread(&mpi_nproc,sizeof(mpi_nproc),1,fp);
        cnt += fread(&mpi_rank,sizeof(mpi_rank),1,fp);
        cnt += fread(&npart_local,sizeof(npart_local),1,fp);
        cnt += fread(&npart_total,sizeof(npart_total),1,fp);

        cnt += fread(&cosm.omegam,sizeof(cosm.omegam),1,fp);
        cnt += fread(&cosm.omegab,sizeof(cosm.omegab),1,fp);
        cnt += fread(&cosm.omegav,sizeof(cosm.omegav),1,fp);
        cnt += fread(&cosm.omeganu,sizeof(cosm.omeganu),1,fp);
        cnt += fread(&cosm.hubble,sizeof(cosm.hubble),1,fp);
        cnt += fread(&anow, sizeof(anow),1,fp);
        cnt += fread(&tnow, sizeof(tnow),1,fp);

        cnt += fread(&lunit,sizeof(lunit),1,fp);
        cnt += fread(&munit,sizeof(munit),1,fp);
        cnt += fread(&tunit,sizeof(tunit),1,fp);

    }

    void update_expansion(COSM::REAL tnow) {
        COSM::REAL om, ov;
        
        om = cosm.omegam;
        ov = cosm.omegav;

        anow = cosm.timetoa(tnow);
        znow = 1.0/anow-1.0;

        hnow = sqrt(1.e0+om*(1.e0/anow-1.e0)+ov*(SQR(anow)-1.e0))/anow;
    }  

    void open_diag_file(char *filename) {
        if (mpi_rank == 0) diag_fp = fopen(filename,"w");
    }

    void output_diag(PS::F64 dtime) {
        static struct timeval prev_tv, now_tv;
        if (mpi_rank == 0) {
            gettimeofday(&now_tv, NULL);
            if (step == 0) {
                fprintf(diag_fp, 
                        "#step   time        dt          redshift    a(t)        wall time\n");
                fprintf(diag_fp,
                        "%5d %11.3e %11.3e %11.3e %11.3e\n",
                        step, tnow, dtime, znow, anow);
            } else {
                float walltime = wallclock_timing(prev_tv, now_tv);
                fprintf(diag_fp,
                        "%5d %11.3e %11.3e %11.3e %11.3e %11.3e\n",
                        step, tnow, dtime, znow, anow, walltime);
            }
            fflush(diag_fp);
            prev_tv = now_tv;
        }
    }

    void close_diag_file(){
        if (mpi_rank == 0) fclose(diag_fp);
    }
};



#endif /* __RUN_PARAM_HPP__ */
