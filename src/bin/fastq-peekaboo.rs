use clap::Parser;
use fastq_peekaboo::{Fastq, FastqRecordRaw};
use std::{
    io::{self, BufWriter, Write},
    path::PathBuf,
};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("File does not exist: {0}")]
    FileDoesNotExist(PathBuf),
    #[error("Provided path is not a file: {0}")]
    NotAFile(PathBuf),
    #[error("{0}")]
    Io(#[from] std::io::Error),
}

#[derive(Parser)]
#[command(version, about = "Fast random-access viewer for FASTQ files")]
struct Args {
    /// Read 1 (or single-ended) FASTQ file path
    #[arg(short = '1', long)]
    r1: PathBuf,

    /// Read 2 FASTQ file path (for paired-end data)
    #[arg(short = '2', long)]
    r2: Option<PathBuf>,

    /// Print the total number of reads
    #[arg(short, long)]
    count: bool,

    /// Fetch a read at this 0-based index
    #[arg(short, long)]
    index: Option<usize>,

    /// Start of range to fetch (0-based, inclusive; requires --end)
    #[arg(short, long, requires = "end")]
    start: Option<usize>,

    /// End of range to fetch (exclusive; requires --start)
    #[arg(short, long, requires = "start")]
    end: Option<usize>,
}

fn main() -> Result<(), Error> {
    let args = Args::parse();

    let fq = if let Some(r2) = args.r2 {
        Fastq::paired(check_path(args.r1)?, check_path(r2)?)?
    } else {
        Fastq::single(check_path(args.r1)?)?
    };

    let mut out = BufWriter::new(io::stdout().lock());

    if args.count {
        if fq.is_paired() {
            writeln!(out, "R1: {} records", fq.read_count_1())?;
            writeln!(out, "R2: {} records", fq.read_count_2())?;
        } else {
            writeln!(out, "{}", fq.read_count_1())?;
        }
        return Ok(());
    }

    if let Some(i) = args.index {
        match fq.raw_read_1(i) {
            Some(rec) => write_fastq_record_raw(&mut out, rec)?,
            None => {
                eprintln!(
                    "Record index {i} out of range (valid range: [0,{}])",
                    fq.read_count_1() - 1
                );
            }
        }
        if let Some(rec) = fq.raw_read_2(i) {
            write_fastq_record_raw(&mut out, rec)?;
        }
    }

    if let (Some(start), Some(end)) = (args.start, args.end) {
        if fq.is_paired() {
            for (r1, r2) in fq.iter_raw_reads_paired(start..end) {
                write_fastq_record_raw(&mut out, r1)?;
                write_fastq_record_raw(&mut out, r2)?;
            }
        } else {
            for rec in fq.iter_raw_reads_1(start..end) {
                write_fastq_record_raw(&mut out, rec)?;
            }
        }
    }

    Ok(())
}

fn write_fastq_record_raw(
    writer: &mut impl Write,
    record: FastqRecordRaw<'_>,
) -> io::Result<()> {
    for line in &record {
        writer.write_all(line.as_bytes())?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

fn check_path(path: PathBuf) -> Result<PathBuf, Error> {
    if !path.exists() {
        return Err(Error::FileDoesNotExist(path));
    }
    if !path.is_file() {
        return Err(Error::NotAFile(path));
    }
    Ok(path)
}
