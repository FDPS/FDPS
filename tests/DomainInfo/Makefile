check:
	make -C collectSampleParticle CCC="$(CCC)" CFLAGS="$(CFLAGS)" PS_PATH="$(PS_PATH)"
	make -C decomposeDomain CCC="$(CCC)" CFLAGS="$(CFLAGS)" PS_PATH="$(PS_PATH)"

clean:
	rm -f *~
	make -C collectSampleParticle clean
	make -C decomposeDomain clean

distclean:
	rm -f *~
	make -C collectSampleParticle distclean
	make -C decomposeDomain distclean


allclean:
	rm -f *~
	make -C collectSampleParticle allclean
	make -C decomposeDomain allclean

