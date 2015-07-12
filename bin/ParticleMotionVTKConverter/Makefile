BIN=ParticleMotionVTKConverter
PF=ParticleMotionVTKConverter.pf
#MAN1=ParticleMotionVTKConverter.1
cxxflags=-g
ldlibs=-lseispp -lgclgrid -lgeocoords -lmwtpp -ltrvltm $(DBLIBS) $(TRLIBS) -lperf -lm
SUBDIR=/contrib
include $(ANTELOPEMAKE)
include $(ANTELOPE)/contrib/include/antelopemake.local
CXXFLAGS += -I$(BOOSTINCLUDE)
LDFLAGS += -L$(BOOSTLIB)

OBJS=ParticleMotionVTKConverter.o ParticleMotionData.o vtp_writer.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
