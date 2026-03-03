/// Errors returned by the YARA mapper.
#[derive(Debug)]
#[non_exhaustive]
pub enum YaraError {
    /// Error opening or loading the YARA FM index.
    IndexOpen(String),
    /// Error building a YARA FM index from a FASTA file.
    IndexBuild(String),
    /// Error during read mapping.
    Mapping(String),
    /// Invalid input data (null bytes, length mismatch, etc.).
    InvalidInput(String),
}

impl std::fmt::Display for YaraError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::IndexOpen(msg) => write!(f, "index open error: {msg}"),
            Self::IndexBuild(msg) => write!(f, "index build error: {msg}"),
            Self::Mapping(msg) => write!(f, "mapping error: {msg}"),
            Self::InvalidInput(msg) => write!(f, "invalid input: {msg}"),
        }
    }
}

impl std::error::Error for YaraError {}
