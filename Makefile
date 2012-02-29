.PHONY: doc
all:
	cd src && make all
doc:
	doxygen doc/Doxyfile
clean:
	cd src && make clean
	rm -rf bin/qsimu
