BIN=test

cflags=-g
cxxflags=-g -I/opt/boost/include
ldflags=-L$(ANTELOPE)/contrib/static -L/opt/boost/lib 
ldlibs=-lmwtpp -lseispp -lgclgrid -lmultiwavelet -lgenloc $(DBLIBS) $(TRLIBS) -lperf -lboost_serialization
SUBDIR=/contrib

include $(ANTELOPEMAKE) 

OBJS=test.o 
$(BIN) : $(OBJS)
	$(RM) $@
	$(CXX) $(CCFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)
