CC = g++
#CFLAGS += -pg
CFLAGS += -O3 -ffast-math
CFLAGS += -Wall
all : plot

plot : plot.o 
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) && mv $@ ../bin/

plot.o : %.o : %.cpp
	$(CC) $(CFLAGS) -c $<

clean :
	rm -rf *.o
