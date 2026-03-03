/// A single CIGAR operation with BAM encoding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarOp {
    /// Operation code: M=0, I=1, D=2, N=3, S=4, H=5, P=6, ==7, X=8.
    pub op: u8,
    /// Length of the operation.
    pub len: u32,
}

impl CigarOp {
    /// Decode from BAM-encoded uint32: `len << 4 | op`.
    #[must_use]
    pub fn from_bam(encoded: u32) -> Self {
        Self { op: (encoded & 0xF) as u8, len: encoded >> 4 }
    }

    /// Operation as a SAM character.
    #[must_use]
    pub fn op_char(&self) -> char {
        match self.op {
            0 => 'M',
            1 => 'I',
            2 => 'D',
            3 => 'N',
            4 => 'S',
            5 => 'H',
            6 => 'P',
            7 => '=',
            8 => 'X',
            _ => '?',
        }
    }
}

impl std::fmt::Display for CigarOp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}", self.len, self.op_char())
    }
}

/// A single alignment record returned by the YARA mapper.
///
/// This is a fully-owned Rust type — all heap data (CIGAR, sequence, quality,
/// XA tag) has been copied from the C++ side and the C++ memory freed.
#[derive(Debug, Clone)]
#[expect(clippy::struct_excessive_bools, reason = "mirrors the C FFI record layout")]
pub struct YaraRecord {
    /// Index of the read pair in the input batch (0-based).
    pub read_pair_index: u32,
    /// Whether this record is for the first read in the pair.
    pub is_read1: bool,

    // Alignment position
    /// Reference contig index.
    pub contig_id: u32,
    /// 0-based leftmost position on the reference.
    pub pos: u32,
    /// Whether the read is mapped to the reverse strand.
    pub is_reverse: bool,
    /// Whether this is a secondary alignment.
    pub is_secondary: bool,
    /// Whether the read is unmapped.
    pub is_unmapped: bool,

    // Alignment quality
    /// Mapping quality.
    pub mapq: u8,
    /// Edit distance (NM tag).
    pub nm: u8,
    /// Number of co-optimal alignments (X0 tag).
    pub x0: u16,
    /// Number of sub-optimal alignments (X1 tag).
    pub x1: u16,

    // Mate info
    /// Mate's reference contig index.
    pub mate_contig_id: u32,
    /// Mate's 0-based position.
    pub mate_pos: u32,
    /// Template length (TLEN).
    pub tlen: i32,
    /// Full SAM flag field.
    pub flag: u16,

    // CIGAR
    /// CIGAR operations (empty for secondaries without `align_secondary`).
    pub cigar: Vec<CigarOp>,

    // Sequence and quality
    /// Read sequence (None for secondary records).
    pub seq: Option<Vec<u8>>,
    /// Base qualities (None for secondary records).
    pub qual: Option<Vec<u8>>,

    /// XA tag string (only when `secondary_mode=Tag`, otherwise None).
    pub xa: Option<String>,
}

impl YaraRecord {
    /// CIGAR string in SAM format (e.g., "50M2I48M").
    #[must_use]
    pub fn cigar_string(&self) -> String {
        use std::fmt::Write;
        let mut s = String::with_capacity(self.cigar.len() * 4);
        for op in &self.cigar {
            write!(s, "{}{}", op.len, op.op_char()).unwrap();
        }
        s
    }
}
