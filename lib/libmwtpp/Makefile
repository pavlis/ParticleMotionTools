LIB=libmwtpp.a
INCLUDE=MWTransform.h PMTimeSeries.h ParticleMotionEllipse.h ParticleMotionError.h
LICENSES=LICENSE_libmwtpp
SUBDIR=/contrib

include $(ANTELOPEMAKE)
OBJS=	# list all objects to put into the library
OBJS=MWTBundle.o MWTransform.o MWTMatrix.o MWTdata.o MWTwaveform.o \
         ParticleMotionEllipse.o ParticleMotionError.o  PMTimeSeries.o
MWTBundle.cc : MWTransform.h
MWTMatrix.cc : MWTransform.h 
MWTdata.cc : MWTransform.h
MWTransform.cc : MWTransform.h
ParticleMotionEllipse.cc : ParticleMotionEllipse.h
ParticleMotionError.cc : ParticleMotionError.h
PMTimeSeries.cc : PMTimeSeries.h ParticleMotionEllipse.h ParticleMotionError.h

$(LIB) : $(OBJS)
	$(RM) $@
	$(CXXAR) $(CXXARFLAGS) $@ $(OBJS)  
	$(RANLIB) $@
