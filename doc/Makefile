all:
	doxygen Doxyfile
	cd Manual && make all
	mv Manual/manual.pdf QSimulator.pdf
clean:
	rm -rf html
	rm -f doxygen.css header.html footer.html
	rm -f QSimulator.pdf
