#include <iostream>
#include <vector>
#include <cassert>

#include "myheader.hpp"

template <class T>
void loadInitialCondition(std::FILE *fp,
                          PS::ParticleSystem<T> & system)
{
    PS::S32 itmp;
    PS::F64 dtmp;

    PS::S32 nloc;

    std::fscanf(fp, "%lf", &dtmp);
    std::fscanf(fp, "%d", &nloc);
    
    system.setNumberOfParticleLocal(nloc);
    
    for(int i = 0; i < nloc; i++) {
        system[i].id     = i;
        std::fscanf(fp, "%d", &itmp);
        std::fscanf(fp, "%lf", &(system[i].mass));
        std::fscanf(fp, "%lf%lf%lf", &(system[i].pos[0]), &(system[i].pos[1]), &(system[i].pos[2]));
        std::fscanf(fp, "%lf%lf%lf", &(system[i].vel[0]), &(system[i].vel[1]), &(system[i].vel[2]));
        system[i].nj = 0;
        system[i].njreal = 0;
    }

    return;
}

template <class T>
void loadInitialConditionFromGreeMFormat(std::FILE *fp,
                                         PS::ParticleSystem<T> & system)
{
    PS::S64 n;
    PS::F32 omg0;

    PS::S32 itmp;
    PS::F32 ftmp;
    PS::F64 dtmp;

    fread(&itmp, sizeof(PS::S32), 1, fp);
    fread(&n,    sizeof(PS::S32), 1, fp);
    fread(&itmp, sizeof(PS::S32), 1, fp);
    fread(&omg0, sizeof(PS::F32), 1, fp);
    fread(&ftmp, sizeof(PS::F32), 1, fp);
    fread(&ftmp, sizeof(PS::F32), 1, fp);
    fread(&ftmp, sizeof(PS::F32), 1, fp);
    fread(&ftmp, sizeof(PS::F32), 1, fp);
    fread(&ftmp, sizeof(PS::F32), 1, fp);
    fread(&ftmp, sizeof(PS::F32), 1, fp);
    fread(&dtmp, sizeof(PS::F64), 1, fp);
    fread(&dtmp, sizeof(PS::F64), 1, fp);
    fread(&dtmp, sizeof(PS::F64), 1, fp);

    PS::S32 n_all;
    PS::S32 n_loc = (PS::S32)n;
    MPI::COMM_WORLD.Allreduce(&n_loc, &n_all, 1, MPI::INT, MPI::SUM);
    PS::F64 mass = 3. * omg0 / (8 * M_PI * (PS::F32)n_all);

    PS::F32 *fcache = new PS::F32[n*3];
    fread(fcache, sizeof(PS::F32), n*3, fp);
    for(PS::S32 i = 0; i < n; i++) {
        system[i].mass   = mass;
        system[i].pos[0] = fcache[3*i+0];
        system[i].pos[1] = fcache[3*i+1];
        system[i].pos[2] = fcache[3*i+2];        
    }
    fread(fcache, sizeof(PS::F32), n*3, fp);
    for(PS::S32 i = 0; i < n; i++) {
        system[i].vel[0] = fcache[3*i+0];
        system[i].vel[1] = fcache[3*i+1];
        system[i].vel[2] = fcache[3*i+2];
    }
    delete(fcache);

    PS::S64 *dcache = new PS::S64[n];
    fread(dcache, sizeof(PS::S64), n, fp);
    for(PS::S32 i = 0; i < n; i++) {
        system[i].id = (PS::S32)dcache[i];
    }
    delete(dcache);

    system.setNumberOfParticleLocal(n);


    return;
}

template <class T>
void loadInitialConditionFromGreeMFormat2(std::FILE *fp,
                                          PS::ParticleSystem<T> & system)
{
    PS::S64 n;
    PS::F32 omg0;

    PS::S32 itmp;
    PS::F32 ftmp;
    PS::F64 dtmp;

    fread(&itmp, sizeof(PS::S32), 1, fp);
    fread(&n,    sizeof(PS::S32), 1, fp);
    fread(&itmp, sizeof(PS::S32), 1, fp);
    fread(&omg0, sizeof(PS::F32), 1, fp);
    fread(&ftmp, sizeof(PS::F32), 1, fp);
    fread(&ftmp, sizeof(PS::F32), 1, fp);
    fread(&ftmp, sizeof(PS::F32), 1, fp);
    fread(&ftmp, sizeof(PS::F32), 1, fp);
    fread(&ftmp, sizeof(PS::F32), 1, fp);
    fread(&ftmp, sizeof(PS::F32), 1, fp);
    fread(&dtmp, sizeof(PS::F64), 1, fp);
    fread(&dtmp, sizeof(PS::F64), 1, fp);
    fread(&dtmp, sizeof(PS::F64), 1, fp);

    PS::S32 n_all;
    PS::S32 n_loc = (PS::S32)n;
    MPI::COMM_WORLD.Allreduce(&n_loc, &n_all, 1, MPI::INT, MPI::SUM);
    PS::F64 mass = 3. * omg0 / (8 * M_PI * (PS::F32)n_all);

    PS::F64 *fcache = new PS::F64[n*3];
    fread(fcache, sizeof(PS::F64), n*3, fp);
    for(PS::S32 i = 0; i < n; i++) {
        system[i].mass   = mass;
        system[i].pos[0] = fcache[3*i+0];
        system[i].pos[1] = fcache[3*i+1];
        system[i].pos[2] = fcache[3*i+2];        
    }
    fread(fcache, sizeof(PS::F64), n*3, fp);
    for(PS::S32 i = 0; i < n; i++) {
        system[i].vel[0] = fcache[3*i+0];
        system[i].vel[1] = fcache[3*i+1];
        system[i].vel[2] = fcache[3*i+2];
    }
    delete(fcache);

    PS::S64 *dcache = new PS::S64[n];
    fread(dcache, sizeof(PS::S64), n, fp);
    for(PS::S32 i = 0; i < n; i++) {
        system[i].id = (PS::S32)dcache[i];
    }
    delete(dcache);

    system.setNumberOfParticleLocal(n);


    return;
}

static PS::S32 nsampave = 100;
static PS::S32 nptclmax = 32768;
static PS::F32 half_len = 1e3;

int main(int argc, char **argv)
{
    PS::Initialize(argc, argv);

    PS::S32 rank = MPI::COMM_WORLD.Get_rank();
    PS::S32 size = MPI::COMM_WORLD.Get_size();
    PS::DomainInfo dinfo;
    PS::ParticleSystem<FP> system, system2;
    PS::TreeForForceLong<Force, EPI, EPJ, PS::MomentMonopolePeriodic, PS::SPJMonopolePeriodic>::WithCutoff treegrav;
    
    dinfo.initialize();
    dinfo.setDomain(2, 2, 2);
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.setPosRootDomain(0., 1.);

    system.initialize();
//    system.setAverageTargetNumberOfSampleParticlePerProcess(nsampave);
    system.createParticle(nptclmax);
    system.setNumberOfParticleLocal(0);

    {
        char filename[64];
        sprintf(filename, "%s.%04d", argv[1], rank);        
        FILE *fp = fopen(filename, "r");
        assert(fp);
        loadInitialConditionFromGreeMFormat2(fp, system);
        fclose(fp);
    }

    PS::F32 theta = 1.0;
//    PS::F32 theta = 0.5;
//    PS::F32 theta = 0.0;
    PS::U32 nleaf = 10;
    PS::U32 ncrit = 300;
    treegrav.initialize(nptclmax, theta, nleaf, ncrit);

/*
    dinfo.collectSampleParticle(system);
    dinfo.decomposeDomain();
*/
    {
        PS::F64ort *dptr = dinfo.getPointerOfPosDomain();
        for(PS::S32 i = 0; i < size; i++) {
            dptr[i].low_[0]  = 0.5 * (PS::S32)(((i >> 2) & 0x1));
            dptr[i].low_[1]  = 0.5 * (PS::S32)(((i >> 1) & 0x1));
            dptr[i].low_[2]  = 0.5 * (PS::S32)(((i >> 0) & 0x1));
            dptr[i].high_[0] = 0.5 * (PS::S32)(((i >> 2) & 0x1) + 1);
            dptr[i].high_[1] = 0.5 * (PS::S32)(((i >> 1) & 0x1) + 1);
            dptr[i].high_[2] = 0.5 * (PS::S32)(((i >> 0) & 0x1) + 1);
        }
    }

    system.exchangeParticle(dinfo);

    // For two kinds of particles
    /*
    {
        PS::S32 nall = system.getNumberOfParticleLocal();
        PS::S32 nhlf = nall / 2;

        system2.initialize();
        system2.setAverageTargetNumberOfSampleParticlePerProcess(nsampave);
        system2.createParticle(nptclmax);        
        PS::S32 nall2 = 0;
        for(PS::S32 i = nhlf; i < nall; i++) {
            system2[nall2] = system[i];
            nall2++;
        }
        system2.setNumberOfParticleLocal(nall2);

        system.setNumberOfParticleLocal(nhlf);

        PS::PM::ParticleMesh pm;
        pm.setDomainInfoParticleMesh(dinfo);
        pm.setParticleParticleMesh(system);
        pm.setParticleParticleMesh(system2, false);
        pm.calcMeshForceOnly();

        nall = system.getNumberOfParticleLocal();
        for(PS::S32 i = 0; i < nall; i++) {
            system[i].apm = pm.getForce(system[i].getPos());
        }
        nall2 = system2.getNumberOfParticleLocal();
        for(PS::S32 i = 0; i < nall2; i++) {
            system2[i].apm = pm.getForce(system2[i].getPos());
        }

        for(PS::S32 i = 0; i < nall2; i++) {
            system[i+nall] = system2[i];
        }
        system.setNumberOfParticleLocal(nall+nall2);
    }
    */

    {
        PS::PM::ParticleMesh pm;
        pm.calcForceAllAndWriteBack(system, dinfo);
    }

#define DEBUG
#ifdef DEBUG
    treegrav.initializeLocalTree(half_len);
    treegrav.setParticleLocalTree(system);
    treegrav.setRootCell(4.0);
    treegrav.mortonSortLocalTreeOnly();
    treegrav.linkCellLocalTreeOnly();    
    treegrav.calcMomentLocalTreeOnly();
    treegrav.exchangeLocalEssentialTree(dinfo);
    treegrav.setLocalEssentialTreeToGlobalTree();
    treegrav.mortonSortGlobalTreeOnly();
    treegrav.linkCellGlobalTreeOnly();
    treegrav.calcMomentGlobalTreeOnly();
    treegrav.makeIPGroup();
    treegrav.calcForceAndWriteBack(calcForceEpEp(), calcForceEpSp(), system);    
#else
    // It doesn't work. Because 'half_len_grav_glb' is too small!
    treegrav.calcForceAllAndWriteBack(calcForceEpEp(), calcForceEpSp(), system, dinfo);
#endif

    {
        char file[64];
        FILE *fp;
        PS::S32 nptcl = system.getNumberOfParticleLocal();
        PS::S32 theta_int = (int)(theta * 10.);

        sprintf(file, "fdps_n32.theta%02d.pos.%04d", theta_int, rank);
        fp = fopen(file, "w");
        for(PS::S32 i = 0; i < nptcl; i++) {
            fprintf(fp, "%8d", system[i].id);
            fprintf(fp, " %+.10e %+.10e %+.10e",
                    system[i].pos.x, system[i].pos.y, system[i].pos.z);
            fprintf(fp, "\n");
        }
        fclose(fp);

        sprintf(file, "fdps_n32.theta%02d.app.%04d", theta_int, rank);
        fp = fopen(file, "w");
        for(PS::S32 i = 0; i < nptcl; i++) {
            fprintf(fp, "%8d", system[i].id);
            fprintf(fp, " %+.10e %+.10e %+.10e",
                    system[i].acc.x, system[i].acc.y, system[i].acc.z);
            fprintf(fp, "\n");
        }
        fclose(fp);

        sprintf(file, "fdps_n32.theta%02d.apm.%04d", theta_int, rank);
        fp = fopen(file, "w");
        for(PS::S32 i = 0; i < nptcl; i++) {
            fprintf(fp, "%8d", system[i].id);
            fprintf(fp, " %+.10e %+.10e %+.10e",
                    system[i].apm.x, system[i].apm.y, system[i].apm.z);
            fprintf(fp, "\n");
        }
        fclose(fp);

        if(rank == 0) {
            printf("tree_theta %02d\n", theta_int);
        }
    }

    /*
    PS::TreeForForceShort<FP, FP, FP>::Gather sph1;
    PS::TreeForForceShort<FP, FP, FP>::Symmetry sph2;
    sph2.copyLocalTree(sph1);
    */

    PS::Finalize();

    return 0;
}
