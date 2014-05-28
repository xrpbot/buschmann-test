EIGEN_INCL = ../3rdparty/eigen-eigen-3.2.0/

cop_ctrl: cop_ctrl.cxx
	g++ -o cop_ctrl -I$(EIGEN_INCL) cop_ctrl.cxx
