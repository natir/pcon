all: test.c ../../target/release/libpcon.a
	gcc -g -Wall test.c ../../target/release/libpcon.a -lpthread -lm -ldl -llzma -lbz2 -o test

../../target/release/libpcon.a:
	cd ../../
	cargo build --release
