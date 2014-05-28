EIGEN_INCL = ../3rdparty/eigen-eigen-3.2.0/

cop_ctrl: cop_ctrl.cxx Hermite.cxx Hermite.h
	g++ -o cop_ctrl -I$(EIGEN_INCL) cop_ctrl.cxx Hermite.cxx

bvp_test: bvp_test.cxx Hermite.cxx Hermite.h
	g++ -o bvp_test -I$(EIGEN_INCL) bvp_test.cxx Hermite.cxx
