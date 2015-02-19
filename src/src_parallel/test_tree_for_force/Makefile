PS_PATH = -I../

#CCC = time mpiFCCpx
#CCC = time FCCpx
#CFLAGS = -Kfast,openmp -NRtrap
#CFLAGS = -Kfast,openmp
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CFLAGS += -Xg
#CFLAGS += -Nfjcex
#CFLAGS += -Krestp=all
#CFLAGS = -Kfast,openmp,restp=all

CCC = mpicxx
#CCC = g++
CFLAGS = -O3 -ffast-math -funroll-loops -Wall
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION
#TOWDIM = PARTICLE_SIMULATOR_TOW_DIMENSION

#all: test_sph.out test_all.out test_integration.out test_sort.out test_make_local_tree.out test_exchange_LET.out 
all: test_nbody_mono.out test_nbody_quad.out test_nbody_mono_geometric_center.out test_nbody_di_geometric_center.out test_nbody_quad_geometric_center.out
#all: test_nbody_quad_geometric_center.out
#all: test_sph.out
#all: test_all.out
#all: test_sph_constN.out
#all: test_all_periodic.out test_all.out

test_all_periodic.out:test_all_periodic.cc
	$(CCC) $(PS_PATH) $(CFLAGS) -o $@ $<

test_sort.out:test_sort.cc 
	$(CCC) $(PS_PATH) $(CFLAGS) -o $@ $<

test_make_local_tree.out:test_make_local_tree.cc 
	$(CCC) $(PS_PATH) $(CFLAGS) -o $@ $<

test_exchange_LET.out:test_exchange_LET.cc
	$(CCC) $(PS_PATH) $(CFLAGS) -o $@ $<

test_integration.out:test_integration.cc
	$(CCC) $(PS_PATH) $(CFLAGS) -o $@ $<

test_all.out:test_all.cc
	$(CCC) $(PS_PATH) $(CFLAGS) -o $@ $<

test_sph.out:test_sph.cc
	$(CCC) $(PS_PATH) $(CFLAGS) -o $@ $<

test_sph_constN.out:test_sph_constN.cc
	$(CCC) $(PS_PATH) $(CFLAGS) -o $@ $<

test_nbody_mono.out:test_nbody.cc
	$(CCC) -DMONOPOLE $(PS_PATH) $(CFLAGS) -o $@ $<

test_nbody_quad.out:test_nbody.cc
	$(CCC) -DQUADRUPOLE $(PS_PATH) $(CFLAGS) -o $@ $<

test_nbody_mono_geometric_center.out:test_nbody.cc
	$(CCC) -DMONOPOLEGEOMETRICCENTER $(PS_PATH) $(CFLAGS) -o $@ $<

test_nbody_di_geometric_center.out:test_nbody.cc
	$(CCC) -DDIPOLEGEOMETRICCENTER $(PS_PATH) $(CFLAGS) -o $@ $<

test_nbody_quad_geometric_center.out:test_nbody.cc
	$(CCC) -DQUADRUPOLEGEOMETRICCENTER $(PS_PATH) $(CFLAGS) -o $@ $<

#clean:
#	rm test_all_periodic.out test_all.out *.s

#clean:
#	rm test_sph_constN.out *.s

#clean:
#	rm test_sph.out *.s

clean:
	rm test_nbody_mono.out test_nbody_quad.out test_nbody_mono_geometric_center.out test_nbody_quad_geometric_center.out *.s

#clean:
#	rm test_all.out *.s

#clean:
#	rm test_integration.out *.s

#clean:
#	rm test_exchange_LET.out *.s
