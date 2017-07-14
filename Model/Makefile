# Must indent with "tab" *not* spaces (damn!)
# NOTE: I changed the LIBS variable from the commented out version to the newer one (because I used a different installation of sundials

TARGET=program_run_trials

#CXX=g++  -O3 -I ../Sundials/include/ -Wno-deprecated
CXX=g++ -fopenmp -std=c++11 -g -O3 -I /usr/local/include/ -Wno-deprecated
#LIBS=-L ../Sundials/lib -l sundials_ida -l sundials_cvode -l 
LIBS=-L /usr/local/lib/ -l sundials_ida -l sundials_cvode -l sundials_nvecserial -l sundials_nvecserial-0 -lm
CC=$(CXX)
SRCS=run_trials_parallel.cpp worm_trials.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o $(TARGET)

check:
	@echo SRCS = $(SRCS)
	@echo OBJS = $(OBJS)

clean:
	rm -rf $(TARGET) $(OBJS)

depend:
	$(CXX) -MM $(SRCS) > .depends.mk

-include .depends.mk
