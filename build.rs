extern crate cbindgen;

use std::env;
use std::process::Command;

fn main() {
    let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();

    cbindgen::Builder::new()
        .with_crate(crate_dir)
        .generate()
        .expect("Unable to generate bindings")
        .write_to_file("dist/pcon.h");

    let output = Command::new("g++")
        .args(&[
            "dist/test_pcon.c",
            "-I",
            "dist/",
            &("target/".to_string() + &env::var("PROFILE").expect("Env error") + "/libpcon.a"),
            "-lpthread",
            "-ldl",
            "-o",
            "dist/test_pcon",
        ])
        .output()
        .expect("failled to build");

    println!("stdout: {}", String::from_utf8_lossy(&output.stdout));
    println!("stderr: {}", String::from_utf8_lossy(&output.stderr));

    if !output.status.success() {
        println!("BUILD FAILLED !!");
    }

    println!("cargo:rerun-if-changed=src/lib.rs");
    println!("cargo:rerun-if-changed=dist/test_pcon.c");
}
