all:	rigid

rigid:
	make -C ./src/rigid/
	cp ./src/rigid/main.out .
flex:
	make -C ./src/flex/
	cp ./src/flex/main.out .
gpu:
	make -C ./src/cuda/
	cp ./src/cuda/gpu.out .
reuse:
	make reuse -C ./src/cuda/
	cp ./src/cuda/reuse.out .
clean:
	make clean -C ./src/rigid/
	make clean -C ./src/flex/
	make clean -C ./src/cuda/
	rm *.out

distclean:
	make dstclean -C ./src/rigid/
	make dstclean -C ./src/flex/
	make dstclean -C ./src/cuda/
	rm *.out
	rm -rf result
