CC = g++
#CFLAGS += -pg
CFLAGS += -O3 -ffast-math
CFLAGS += -Wall
INC_DIR += /opt/include/
LIB_DIR += /opt/lib
LIBS += -lcvm_em64t -loptions -lmathexpr
all : qsimu

qsimu : main.o 
	$(CC) $(CFLAGS) -I$(INC_DIR) -L$(LIB_DIR) $^ -o $@ $(LIBS) && mv $@ ../bin/

main.o : %.o : %.cpp
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $<

clean :
	rm -rf *.o
