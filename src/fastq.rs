use memmap2::{Advice, Mmap};
use std::{
    fs::File,
    ops::Range,
    path::{Path, PathBuf},
};

pub type FastqRecordRaw<'a> = [&'a str; 4];

#[derive(Debug)]
pub struct Fastq {
    r1: Option<MappedFastqFile>,
    r2: Option<MappedFastqFile>,
    files: FastqFilePaths,
}

#[derive(Debug, Clone)]
pub enum FastqFilePaths {
    Single { f: PathBuf },
    Paired { f1: PathBuf, f2: PathBuf },
}

/// A single memory-mapped FASTQ file with a byte-offset index.
#[derive(Debug)]
struct MappedFastqFile {
    mmap: Mmap,
    offsets: Vec<usize>,
}

impl Fastq {
    pub fn single(f: PathBuf) -> std::io::Result<Self> {
        Ok(Self {
            r1: Some(MappedFastqFile::open(&f)?),
            r2: None,
            files: FastqFilePaths::Single { f },
        })
    }

    pub fn paired(f1: PathBuf, f2: PathBuf) -> std::io::Result<Self> {
        Ok(Self {
            r1: Some(MappedFastqFile::open(&f1)?),
            r2: Some(MappedFastqFile::open(&f2)?),
            files: FastqFilePaths::Paired { f1, f2 },
        })
    }

    pub fn files(&self) -> FastqFilePaths {
        self.files.clone()
    }

    pub fn is_paired(&self) -> bool {
        matches!(self.files, FastqFilePaths::Paired { .. })
    }

    pub fn read_count_1(&self) -> usize {
        self.r1.as_ref().map_or(0, |m| m.read_count())
    }

    pub fn read_count_2(&self) -> usize {
        self.r2.as_ref().map_or(0, |m| m.read_count())
    }

    pub fn raw_read_1(&self, i: usize) -> Option<FastqRecordRaw<'_>> {
        self.r1.as_ref()?.get(i)
    }

    pub fn raw_read_2(&self, i: usize) -> Option<FastqRecordRaw<'_>> {
        self.r2.as_ref()?.get(i)
    }

    pub fn raw_read_pair(
        &self,
        i: usize,
    ) -> Option<(FastqRecordRaw<'_>, FastqRecordRaw<'_>)> {
        Some((self.raw_read_1(i)?, self.raw_read_2(i)?))
    }

    pub fn iter_raw_reads_1(
        &self,
        range: Range<usize>,
    ) -> impl Iterator<Item = FastqRecordRaw<'_>> {
        (range.start..range.end).filter_map(move |i| self.raw_read_1(i))
    }

    pub fn iter_raw_reads_2(
        &self,
        range: Range<usize>,
    ) -> impl Iterator<Item = FastqRecordRaw<'_>> {
        (range.start..range.end).filter_map(move |i| self.raw_read_2(i))
    }

    pub fn iter_raw_reads_paired(
        &self,
        range: Range<usize>,
    ) -> impl Iterator<Item = (FastqRecordRaw<'_>, FastqRecordRaw<'_>)> {
        (range.start..range.end).filter_map(move |i| self.raw_read_pair(i))
    }
}

impl MappedFastqFile {
    /// Return the four lines of a read record `i` as borrowed `&str` slices into the mmap.
    fn get(&self, i: usize) -> Option<FastqRecordRaw<'_>> {
        let &start = self.offsets.get(i)?;
        let end = self.offsets.get(i + 1).copied().unwrap_or(self.byte_count());
        // Skipping validity checks for speed. Think how likely this is to fail.
        let text =
            unsafe { std::str::from_utf8_unchecked(&self.mmap[start..end]) };
        let mut lines = text.lines();
        Some([lines.next()?, lines.next()?, lines.next()?, lines.next()?])
    }

    // Memory-map the file and build the record-offset index.
    // This is the only O(n) operation; every access after this is O(1).
    fn open(path: &Path) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };

        #[cfg(unix)]
        let _ = mmap.advise(Advice::Sequential);

        let offsets = Self::build_index(&mmap);

        #[cfg(unix)]
        let _ = mmap.advise(Advice::Random);

        Ok(Self { mmap, offsets })
    }

    fn build_index(bytes: &[u8]) -> Vec<usize> {
        if bytes.is_empty() {
            return Vec::new();
        }

        let mut offsets = Vec::with_capacity(bytes.len() / 400);

        offsets.push(0);
        let mut line_count: u64 = 0;
        for pos in memchr::memchr_iter(b'\n', bytes) {
            line_count += 1;
            if line_count & 3 == 0 {
                let next = pos + 1;
                if next < bytes.len() {
                    offsets.push(next);
                }
            }
        }

        offsets
    }

    #[inline]
    fn read_count(&self) -> usize {
        self.offsets.len()
    }

    #[inline]
    fn byte_count(&self) -> usize {
        self.mmap.len()
    }
}
