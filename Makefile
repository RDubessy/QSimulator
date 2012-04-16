.PHONY: doc tools
all:
	cd src && make all
tools:
	cd tools && make all
doc:
	doxygen doc/Doxyfile
clean:
	cd src && make clean
	cd tools && make clean
	rm -rf bin/qsimu
