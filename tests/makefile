copy:
	cp -r ../../../tests/test1 ./
	cp -r ../../../tests/test2 ./
	cp -r ../../../tests/test3 ./
	cp -r ../../../tests/test4 ./

testbase:
	make -C 01_steady_flow_123d test	

testall: 
	make -C 01_steady_flow_123d test
	make -C 02_transport_12d test
	make -C 03_transport_small_12d test
	make -C 04_density_driven_12d test

clean:
	rm -rf ./test1
	rm -rf ./test2
	rm -rf ./test3
	rm -rf ./test4
	
%.tst:
	@BASE=$*;\
	dir="$${BASE}*";\
	if [ -d $${dir} ];\
	then make -C $${dir} test;\
	else echo "missing test directory $${dir}";\
	fi