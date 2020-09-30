all: test.c ../../target/debug/libpcon.a
	gcc -g -Wall test.c ../../target/debug/libpcon.a -lpthread -lm -ldl -llzma -lbz2 -o test

../../target/debug/libpcon.a:
	cargo build --manifest-path ../../Cargo.toml
