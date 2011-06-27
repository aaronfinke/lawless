noinst_LIBRARIES = libbfgs.a

libbfgs_a_SOURCES = RefineBase.cc RefineBase2.cc Minimizer.cc PhaserError.cc Errors.cc

noinst_HEADERS = Errors.h        PhaserError.h   RefineBase2.h   maths.h         symm_eigen.h \
Minimizer.h     ProtocolBase.h  fmat.h          phaser_types.h  vec.h \
Output.h        RefineBase.h    jiffy.h         subscript.h


# Dependencies:
#  CCTBX
#  TNT

INCLUDES = $(PTHREAD_CFLAGS) $(CCTBX_CXXFLAGS) 

