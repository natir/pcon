//! Command Line Interface declaration of project pcon

/* std use */
use std::io::Read as _;

/* crate use */

/* project use */
use crate::error;

/// Prompt COuNter, a short kmer counter.
#[derive(clap::Parser, std::fmt::Debug)]
#[clap(
    name = "pcon",
    version = "0.1",
    author = "Pierre Marijon <pierre@marijon.fr>"
)]
#[clap(propagate_version = true)]
pub struct Command {
    /// SubCommand
    #[clap(subcommand)]
    pub subcommand: SubCommand,

    #[cfg(feature = "parallel")]
    /// Number of theard use 0 use all avaible core, default value 0
    #[clap(short = 't', long = "threads")]
    threads: Option<usize>,

    /// Silence all output
    #[clap(short = 'q', long = "quiet")]
    quiet: bool,

    /// Verbose mode (-v, -vv, -vvv, etc)
    #[clap(short = 'v', long = "verbosity", parse(from_occurrences))]
    verbosity: usize,

    /// Timestamp (sec, ms, ns, none)
    #[clap(short = 'T', long = "timestamp")]
    ts: Option<stderrlog::Timestamp>,
}

impl Command {
    /// Get number of thread
    #[cfg(feature = "parallel")]
    pub fn threads(&self) -> usize {
        self.threads.unwrap_or(0)
    }

    /// Get verbosity level
    pub fn verbosity(&self) -> usize {
        self.verbosity
    }

    /// Get quiet
    pub fn quiet(&self) -> bool {
        self.quiet
    }

    /// Get timestamp granularity
    pub fn timestamp(&self) -> stderrlog::Timestamp {
        self.ts.unwrap_or(stderrlog::Timestamp::Off)
    }
}

/// Enumeration of subcommand
#[derive(clap::Subcommand, std::fmt::Debug)]
pub enum SubCommand {
    /// Perform count of kmer
    Count(Count),

    /// Convert pcon native output in other type
    Dump(Dump),
}

/// Choose dump type
#[derive(Copy, Clone, Eq, Debug, PartialEq, PartialOrd, Ord, clap::ValueEnum)]
pub enum DumpType {
    /// Output in bin mode
    Pcon,

    /// Output in csv mode
    Csv,

    /// Output in solid mode
    Solid,

    /// Output in spectrum mode
    Spectrum,
}

impl std::str::FromStr for DumpType {
    type Err = error::Error;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "pcon" => Ok(DumpType::Pcon),
            "csv" => Ok(DumpType::Csv),
            "solid" => Ok(DumpType::Solid),
            "spectrum" => Ok(DumpType::Spectrum),
            _ => Err(error::Error::DumpTypeFromStr(s.to_string())),
        }
    }
}

/// SubCommand Count
#[derive(clap::Args, std::fmt::Debug)]
pub struct Count {
    /// Size of kmer
    #[clap(short = 'k', long = "kmer-size")]
    kmer_size: u8,

    /// Path to inputs, default read stdin
    #[clap(short = 'i', long = "inputs")]
    inputs: Option<Vec<std::path::PathBuf>>,

    /// Path where count are store, default write in stdout
    #[clap(short = 'o', long = "output")]
    output: Option<std::path::PathBuf>,

    /// Minimal abundance, default value 0
    #[clap(short = 'a', long = "abundance")]
    abundance: Option<u8>,

    /// Dump type, default bin
    #[clap(short = 'd', long = "dump")]
    dump: Option<DumpType>,

    /// Number of sequence record load in buffer, default 8192
    #[clap(short = 'b', long = "record_buffer")]
    record_buffer: Option<u64>,
}

impl Count {
    /// Get size of kmer
    pub fn kmer_size(&self) -> u8 {
        self.kmer_size - (!(self.kmer_size & 0b1) & 0b1)
    }

    /// Get inputs
    pub fn inputs(&self) -> error::Result<Box<dyn std::io::BufRead>> {
        match &self.inputs {
            None => Ok(Box::new(std::io::stdin().lock())),
            Some(paths) => {
                let mut handle: Box<dyn std::io::Read> = Box::new(std::io::Cursor::new(vec![]));

                for path in paths {
                    let (file, _compression) =
                        niffler::get_reader(Box::new(std::fs::File::open(path)?))?;
                    handle = Box::new(handle.chain(file));
                }

                Ok(Box::new(std::io::BufReader::new(handle)))
            }
        }
    }

    /// Get output
    pub fn output(&self) -> error::Result<Box<dyn std::io::Write>> {
        match &self.output {
            None => Ok(Box::new(std::io::BufWriter::new(std::io::stdout().lock()))),
            Some(path) => {
                let handle: Box<dyn std::io::Write> =
                    Box::new(std::fs::File::create(path).map(std::io::BufWriter::new)?);

                Ok(handle)
            }
        }
    }

    /// Get dump type
    pub fn dump(&self) -> DumpType {
        match self.dump {
            Some(a) => a,
            None => DumpType::Pcon,
        }
    }

    /// Get abundance
    pub fn abundance(&self) -> u8 {
        self.abundance.unwrap_or(0)
    }

    /// Get record_buffer
    pub fn record_buffer(&self) -> u64 {
        self.record_buffer.unwrap_or(8192)
    }
}

/// SubCommand Dump
#[derive(clap::Args, std::fmt::Debug)]
pub struct Dump {
    /// Path to inputs, default read stdin
    #[clap(short = 'i', long = "inputs")]
    input: Option<std::path::PathBuf>,

    /// Path where count are store, default write in stdout
    #[clap(short = 'o', long = "output")]
    output: Option<std::path::PathBuf>,

    /// Dump type
    #[clap(value_enum)]
    dump: DumpType,

    /// Minimal abundance, default value 0
    #[clap(short = 'a', long = "abundance")]
    abundance: u8,
}

impl Dump {
    /// Get inputs
    pub fn input(&self) -> error::Result<Box<dyn std::io::BufRead>> {
        match &self.input {
            None => Ok(Box::new(std::io::stdin().lock())),
            Some(path) => {
                let handle: Box<dyn std::io::Read> = Box::new(std::fs::File::open(path)?);

                Ok(Box::new(std::io::BufReader::new(handle)))
            }
        }
    }

    /// Get output
    pub fn output(&self) -> error::Result<Box<dyn std::io::Write>> {
        match &self.output {
            None => Ok(Box::new(std::io::BufWriter::new(std::io::stdout().lock()))),
            Some(path) => {
                let handle: Box<dyn std::io::Write> =
                    Box::new(std::fs::File::create(path).map(std::io::BufWriter::new)?);

                Ok(handle)
            }
        }
    }

    /// Get dump type
    pub fn dump(&self) -> DumpType {
        self.dump
    }

    /// Get abundance
    pub fn abundance(&self) -> u8 {
        self.abundance
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::io::Write as _;
    use std::str::FromStr as _;

    #[test]
    fn dumptype_conversion() -> error::Result<()> {
        assert_eq!(DumpType::from_str("pcon")?, DumpType::Pcon);
        assert_eq!(DumpType::from_str("Pcon")?, DumpType::Pcon);
        assert_eq!(DumpType::from_str("PCON")?, DumpType::Pcon);

        assert_eq!(DumpType::from_str("csv")?, DumpType::Csv);
        assert_eq!(DumpType::from_str("Csv")?, DumpType::Csv);
        assert_eq!(DumpType::from_str("CSV")?, DumpType::Csv);

        assert_eq!(DumpType::from_str("solid")?, DumpType::Solid);
        assert_eq!(DumpType::from_str("Solid")?, DumpType::Solid);
        assert_eq!(DumpType::from_str("SOLID")?, DumpType::Solid);

        assert_eq!(DumpType::from_str("spectrum")?, DumpType::Spectrum);
        assert_eq!(DumpType::from_str("Spectrum")?, DumpType::Spectrum);
        assert_eq!(DumpType::from_str("SPECTRUM")?, DumpType::Spectrum);

        assert!(DumpType::from_str("").is_err());

        Ok(())
    }

    #[cfg(not(feature = "parallel"))]
    #[test]
    fn basic() {
        let subcmd = Count {
            inputs: None,
            output: None,
            kmer_size: 32,
            abundance: Some(0),
            dump: None,
            record_buffer: None,
        };

        let cmd = Command {
            verbosity: 3,
            quiet: false,
            ts: None,
            subcommand: SubCommand::Count(subcmd),
        };

        assert_eq!(cmd.verbosity(), 3);
        assert!(!cmd.quiet());
        assert!(matches!(cmd.timestamp(), stderrlog::Timestamp::Off));
    }

    #[cfg(feature = "parallel")]
    #[test]
    fn basic_parallel() {
        let subcmd = Count {
            inputs: None,
            output: None,
            kmer_size: 32,
            abundance: None,
            dump: None,
            record_buffer: None,
        };

        let cmd = Command {
            verbosity: 3,
            quiet: false,
            ts: None,
            subcommand: SubCommand::Count(subcmd),
            threads: Some(8),
        };

        assert_eq!(cmd.verbosity(), 3);
        assert!(!cmd.quiet());
        assert!(matches!(cmd.timestamp(), stderrlog::Timestamp::Off));
        assert_eq!(cmd.threads(), 8);
    }

    #[test]
    fn count() -> error::Result<()> {
        let mut input1 = tempfile::NamedTempFile::new()?;
        input1.write_all(b">test\n")?;
        let mut input2 = tempfile::NamedTempFile::new()?;
        input2.write_all(b"TACG\n")?;

        let count = Count {
            inputs: Some(vec![
                input1.path().to_path_buf(),
                input2.path().to_path_buf(),
            ]),
            output: None,
            kmer_size: 32,
            abundance: Some(2),
            dump: Some(DumpType::Solid),
            record_buffer: Some(512),
        };

        let mut content = Vec::new();
        count.inputs()?.read_to_end(&mut content)?;
        assert_eq!(content, b">test\nTACG\n");

        assert_eq!(count.kmer_size(), 31);
        assert_eq!(count.abundance(), 2);
        assert_eq!(count.dump(), DumpType::Solid);
        assert_eq!(count.record_buffer(), 512);

        Ok(())
    }

    #[test]
    fn dump() -> error::Result<()> {
        let mut input1 = tempfile::NamedTempFile::new()?;
        input1.write_all(b">test\nATCG\n")?;

        let dump = Dump {
            input: Some(input1.path().to_path_buf()),
            output: None,
            abundance: 2,
            dump: DumpType::Solid,
        };

        let mut content = Vec::new();
        dump.input()?.read_to_end(&mut content)?;
        assert_eq!(content, b">test\nATCG\n");

        assert_eq!(dump.abundance(), 2);
        assert_eq!(dump.dump(), DumpType::Solid);

        Ok(())
    }
}
