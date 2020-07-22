if [ ! -d build ]
then
	mkdir build
fi
cd build
cmake ..
make
if [ ! -z $TESTS ]
then
	ctest -C Debug
fi
make install
