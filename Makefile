.PHONY: doc tools
all:
	cd src && make all
	make tools
install:
	mv bin/qsimu /usr/local/bin/qsimu
tools:
	cd tools && make all
doc:
	cd doc && make all
clean:
	cd src && make clean
	cd tools && make clean
	cd doc && make clean
	rm -rf bin/*
