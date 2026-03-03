/// Sensitivity level for YARA's seed-and-extend pipeline.
///
/// - `Low`: Only two rounds of seeding (0-error and 1-error seeds).
/// - `High`: Three rounds of seeding (adds 2-error seeds). Default.
/// - `Full`: Same as `High` but uses edit-distance seeds instead of Hamming.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum Sensitivity {
    Low = 0,
    High = 1,
    Full = 2,
}

/// How secondary alignments are reported.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum SecondaryMode {
    /// Encode secondary alignments in the XA tag of the primary record.
    Tag = 0,
    /// Emit secondary alignments as separate records with the SECONDARY flag.
    Record = 1,
    /// Do not report secondary alignments.
    Omit = 2,
}

/// Configuration options for the YARA mapper.
///
/// Use [`MapperOptions::default()`] for sensible defaults matching YARA's
/// command-line defaults, then override individual fields as needed.
#[derive(Debug, Clone)]
pub struct MapperOptions {
    /// Maximum error rate as a fraction of read length (default 0.05 = 5%).
    pub error_rate: f32,
    /// Strata rate — fraction of additional errors beyond the best alignment
    /// to report (default 0.00).
    pub strata_rate: f32,
    /// Fixed strata count (alternative to `strata_rate`). Set to `None` to
    /// use `strata_rate` instead.
    pub strata_count: Option<u32>,
    /// Sensitivity level (default `High`).
    pub sensitivity: Sensitivity,
    /// Number of threads for OpenMP parallelism (default: available cores).
    pub threads: usize,
    /// How to report secondary alignments (default `Tag`).
    pub secondary_mode: SecondaryMode,
    /// Whether to align secondary matches (compute CIGAR). Only meaningful
    /// when `secondary_mode` is `Tag` or `Record`.
    pub align_secondary: bool,
    /// Whether to verify mate-pair matches (paired-end rescue). Default true.
    pub verify_matches: bool,
    /// Verbosity level (0 = silent, 1 = stats, 2 = detailed).
    pub verbose: u32,
    /// Expected library (insert) length. 0 = auto-estimate from the data.
    pub library_length: u32,
    /// Expected library standard deviation. 0 = auto-estimate.
    pub library_dev: u32,
}

impl Default for MapperOptions {
    fn default() -> Self {
        Self {
            error_rate: 0.05,
            strata_rate: 0.00,
            strata_count: None,
            sensitivity: Sensitivity::High,
            threads: std::thread::available_parallelism().map(std::num::NonZero::get).unwrap_or(1),
            secondary_mode: SecondaryMode::Tag,
            align_secondary: false,
            verify_matches: true,
            verbose: 0,
            library_length: 0,
            library_dev: 0,
        }
    }
}

impl MapperOptions {
    #[expect(
        clippy::cast_possible_wrap,
        clippy::cast_possible_truncation,
        reason = "FFI boundary: all values fit in i32 for the C API"
    )]
    pub(crate) fn to_ffi(&self) -> yara_mapper_sys::YaraMapperOptions {
        yara_mapper_sys::YaraMapperOptions {
            error_rate: self.error_rate,
            strata_rate: self.strata_rate,
            strata_count: self.strata_count.map_or(-1, |c| c as i32),
            sensitivity: self.sensitivity as i32,
            threads: self.threads as i32,
            secondary_mode: self.secondary_mode as i32,
            align_secondary: i32::from(self.align_secondary),
            verify_matches: i32::from(self.verify_matches),
            verbose: self.verbose as i32,
            library_length: self.library_length as i32,
            library_dev: self.library_dev as i32,
        }
    }
}
