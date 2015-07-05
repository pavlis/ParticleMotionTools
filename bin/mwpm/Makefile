BIN=test

cflags=-g
cxxflags=-g
ldflags=-L$(ANTELOPE)/contrib/static
ldlibs=-lmwtpp -lseispp -lgclgrid -lmultiwavelet -lgenloc $(DBLIBS) $(TRLIBS) -lperf
SUBDIR=/contrib

include $(ANTELOPEMAKE) 

OBJS=PMTseriesToVTK.o vtp_writer.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
