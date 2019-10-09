all: io_raw_cc copy_file

io_raw_cc:
	g++ io_raw_cc.cpp -o io_raw_cc

copy_file:
	mv -f io_raw_cc cc_test

