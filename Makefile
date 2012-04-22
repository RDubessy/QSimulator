.PHONY: doc tools
all:
	cd src && make all
	make tools
	make doc
tools:
	cd tools && make all
doc:
	doxygen doc/Doxyfile
clean:
	cd src && make clean
	cd tools && make clean
	rm -rf bin/*
