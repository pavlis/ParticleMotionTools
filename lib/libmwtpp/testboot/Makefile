BIN=testboot
cxxflags=-g -I/opt/boost/include
ldlibs=-L.. -L/opt/boost/lib -lseispp -lgclgrid -lmwtpp  $(DBLIBS) $(TRLIBS) -lperf -lboost_serialization -lseispp 

SUBDIR=/contrib

include $(ANTELOPEMAKE)  	
include $(ANTELOPEMAKELOCAL)

OBJS=testboot.o
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
