all:  build install
      
        
build: build/CMakeCache.txt
	make -C build tpeg

        
# creates out of source build directory
build/CMakeCache.txt:
	if [ ! -d build ]; then mkdir build; fi
	cd build; cmake ..

        
install:
	cp build/bin/tpeg bin/tpeg
	
	

clean-all:
	rm -rf build
	rm -f bin/tpeg