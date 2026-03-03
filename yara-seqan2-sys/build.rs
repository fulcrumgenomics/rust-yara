use std::path::PathBuf;

fn main() {
    let manifest_dir =
        PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR not set"));

    let include_dir = manifest_dir.join("vendor").join("include");
    let yara_app_dir = manifest_dir.join("vendor").join("apps").join("yara");

    assert!(
        include_dir.is_dir(),
        "vendored SeqAn2 include/ not found at {}",
        include_dir.display()
    );
    assert!(yara_app_dir.is_dir(), "vendored YARA app/ not found at {}", yara_app_dir.display());

    println!("cargo::metadata=include={}", include_dir.display());
    println!("cargo::metadata=yara_app={}", yara_app_dir.display());
}
