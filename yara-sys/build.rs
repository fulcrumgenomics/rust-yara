use std::env;
use std::path::PathBuf;

fn main() {
    // Use SEQAN_DIR as an override for local development; otherwise read paths
    // from the vendored yara-seqan2-sys crate via cargo's DEP_SEQAN2_* env vars.
    let (seqan_include, yara_app_dir) = if let Ok(seqan_dir) = env::var("SEQAN_DIR") {
        (PathBuf::from(&seqan_dir).join("include"), PathBuf::from(&seqan_dir).join("apps/yara"))
    } else {
        let include = env::var("DEP_SEQAN2_INCLUDE")
            .expect("DEP_SEQAN2_INCLUDE not set (is yara-seqan2-sys a dependency?)");
        let yara_app = env::var("DEP_SEQAN2_YARA_APP")
            .expect("DEP_SEQAN2_YARA_APP not set (is yara-seqan2-sys a dependency?)");
        (PathBuf::from(include), PathBuf::from(yara_app))
    };

    // Detect OpenMP on macOS via Homebrew libomp.
    let (omp_include, omp_lib) = detect_openmp();

    let mut build = cc::Build::new();
    build
        .cpp(true)
        .file("cpp/yara_shim.cpp")
        .file("cpp/yara_indexer_shim.cpp")
        .include(seqan_include)
        .include(&yara_app_dir) // YARA app headers
        .include("cpp") // shim headers
        .std("c++17")
        .opt_level_str("2")
        .define("SEQAN_HAS_ZLIB", "1")
        .define("YARA_LARGE_CONTIGS", "1")
        // Suppress common SeqAn2 warnings that are harmless.
        .flag("-Wno-unused-parameter")
        .flag("-Wno-sign-compare")
        .flag("-Wno-unused-variable")
        .flag("-Wno-deprecated-declarations");

    // OpenMP support
    if let Some(ref inc) = omp_include {
        build.flag(format!("-I{}", inc.display()));
        build.flag("-Xpreprocessor");
        build.flag("-fopenmp");
    } else {
        build.flag("-fopenmp");
    }

    build.compile("yara_shim");

    // Link dependencies
    println!("cargo:rustc-link-lib=z");

    if let Some(ref lib) = omp_lib {
        println!("cargo:rustc-link-search=native={}", lib.display());
        println!("cargo:rustc-link-lib=omp");
    } else {
        println!("cargo:rustc-link-lib=gomp");
    }

    // Link C++ standard library
    let target = env::var("TARGET").unwrap_or_default();
    if target.contains("apple") {
        println!("cargo:rustc-link-lib=c++");
    } else {
        println!("cargo:rustc-link-lib=stdc++");
    }

    // Rebuild if shim source changes
    println!("cargo:rerun-if-changed=cpp/yara_shim.h");
    println!("cargo:rerun-if-changed=cpp/yara_shim.cpp");
    println!("cargo:rerun-if-changed=cpp/yara_indexer_shim.h");
    println!("cargo:rerun-if-changed=cpp/yara_indexer_shim.cpp");
    println!("cargo:rerun-if-env-changed=SEQAN_DIR");
}

/// Detect OpenMP (libomp) on macOS via Homebrew, or fall back to system OpenMP.
fn detect_openmp() -> (Option<PathBuf>, Option<PathBuf>) {
    // Try Homebrew libomp
    let homebrew_prefix =
        env::var("HOMEBREW_PREFIX").unwrap_or_else(|_| "/opt/homebrew".to_string());
    let omp_prefix = PathBuf::from(&homebrew_prefix).join("opt").join("libomp");

    let omp_include = omp_prefix.join("include");
    let omp_lib = omp_prefix.join("lib");

    if omp_include.join("omp.h").exists() && omp_lib.join("libomp.dylib").exists() {
        return (Some(omp_include), Some(omp_lib));
    }

    // Fall back: assume system-provided OpenMP (Linux with libgomp)
    (None, None)
}
