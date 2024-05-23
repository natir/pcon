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
    #[clap(short = 'v', long = "verbosity", action = clap::ArgAction::Count)]
    verbosity: u8,

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
        self.verbosity as usize
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

    #[cfg(not(feature = "parallel"))]
    /// Perform count of large kmer if associate minimizer is present more than abundance
    MiniCount(MiniCount),

    /// Convert pcon native output in other format
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
}

/// Choose input format
#[derive(Copy, Clone, Eq, Debug, PartialEq, PartialOrd, Ord, clap::ValueEnum)]
pub enum Format {
    /// Input in format fasta
    Fasta,

    #[cfg(feature = "fastq")]
    /// Input in format fastq
    Fastq,
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

    /// Format of input, default fasta
    #[clap(short = 'f', long = "formats")]
    format: Option<Format>,

    /// Path where count are store, default write in stdout
    #[clap(short = 'p', long = "pcon")]
    pcon: Option<Vec<std::path::PathBuf>>,

    /// Path where count are store
    #[clap(short = 'c', long = "csv")]
    csv: Option<Vec<std::path::PathBuf>>,

    /// Path where count are store
    #[clap(short = 's', long = "solid")]
    solid: Option<Vec<std::path::PathBuf>>,

    /// Minimal abundance, default value 0
    #[clap(short = 'a', long = "abundance")]
    abundance: Option<crate::CountTypeNoAtomic>,

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

    /// Get format inputs
    pub fn format(&self) -> Format {
        self.format.unwrap_or(Format::Fasta)
    }

    /// Get output
    pub fn outputs(
        &self,
    ) -> Vec<(
        DumpType,
        error::Result<Box<dyn std::io::Write + std::marker::Send>>,
    )> {
        let mut outputs: Vec<(
            DumpType,
            error::Result<Box<dyn std::io::Write + std::marker::Send>>,
        )> = vec![];

        match &self.csv {
            None => (),
            Some(paths) => {
                for path in paths {
                    outputs.push((DumpType::Csv, create(path)));
                }
            }
        }

        match &self.solid {
            None => (),
            Some(paths) => {
                for path in paths {
                    outputs.push((DumpType::Solid, create(path)));
                }
            }
        }

        match &self.pcon {
            None => {
                if outputs.is_empty() {
                    outputs.push((
                        DumpType::Pcon,
                        Ok(Box::new(std::io::BufWriter::new(std::io::stdout()))),
                    ))
                }
            }
            Some(paths) => {
                for path in paths {
                    outputs.push((DumpType::Pcon, create(path)));
                }
            }
        }

        outputs
    }

    /// Get abundance
    pub fn abundance(&self) -> crate::CountTypeNoAtomic {
        self.abundance.unwrap_or(0)
    }

    /// Get record_buffer
    pub fn record_buffer(&self) -> u64 {
        self.record_buffer.unwrap_or(8192)
    }
}

/// SubCommand MiniCount
#[derive(clap::Args, std::fmt::Debug)]
pub struct MiniCount {
    /// Size of kmer
    #[clap(short = 'k', long = "kmer-size")]
    kmer_size: u64,

    /// Size of minimizer
    #[clap(short = 'm', long = "minimizer-size")]
    minimizer_size: u8,

    /// Path to inputs, default read stdin
    #[clap(short = 'i', long = "inputs")]
    inputs: Option<Vec<std::path::PathBuf>>,

    /// Format of input, default fasta
    #[clap(short = 'f', long = "formats")]
    format: Option<Format>,

    /// Path where count are store
    #[clap(short = 'c', long = "csv")]
    csv: Option<Vec<std::path::PathBuf>>,

    /// Minimal abundance, default value 0
    #[clap(short = 'a', long = "abundance")]
    abundance: Option<crate::CountTypeNoAtomic>,

    /// Minimal minimizer abundance, default value 2
    #[clap(short = 'A', long = "mini-abundance")]
    mini_abundance: Option<crate::CountTypeNoAtomic>,

    /// Number of sequence record load in buffer, default 8192
    #[clap(short = 'b', long = "record_buffer")]
    record_buffer: Option<u64>,
}

impl MiniCount {
    /// Get size of kmer
    pub fn kmer_size(&self) -> u64 {
        self.kmer_size
    }

    /// Get size of minimizer
    pub fn minimizer_size(&self) -> u8 {
        self.minimizer_size - (!(self.minimizer_size & 0b1) & 0b1)
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

    /// Get format inputs
    pub fn format(&self) -> Format {
        self.format.unwrap_or(Format::Fasta)
    }

    /// Get output
    pub fn outputs(
        &self,
    ) -> Vec<(
        DumpType,
        error::Result<Box<dyn std::io::Write + std::marker::Send>>,
    )> {
        let mut outputs: Vec<(
            DumpType,
            error::Result<Box<dyn std::io::Write + std::marker::Send>>,
        )> = vec![];

        match &self.csv {
            None => outputs.push((
                DumpType::Csv,
                Ok(Box::new(std::io::BufWriter::new(std::io::stdout()))),
            )),
            Some(paths) => {
                for path in paths {
                    outputs.push((DumpType::Csv, create(path)));
                }
            }
        }

        outputs
    }

    /// Get abundance
    pub fn abundance(&self) -> crate::CountTypeNoAtomic {
        self.abundance.unwrap_or(0)
    }

    /// Get abundance
    pub fn mini_abundance(&self) -> crate::CountTypeNoAtomic {
        self.mini_abundance.unwrap_or(2)
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
    #[clap(short = 'c', long = "csv")]
    csv: Option<Vec<std::path::PathBuf>>,

    /// Path where count are store
    #[clap(short = 'p', long = "pcon")]
    pcon: Option<Vec<std::path::PathBuf>>,

    /// Path where count are store
    #[clap(short = 's', long = "solid")]
    solid: Option<Vec<std::path::PathBuf>>,

    /// Minimal abundance, default value 0
    #[clap(short = 'a', long = "abundance")]
    abundance: crate::CountTypeNoAtomic,
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
    pub fn outputs(
        &self,
    ) -> Vec<(
        DumpType,
        error::Result<Box<dyn std::io::Write + std::marker::Send>>,
    )> {
        let mut outputs: Vec<(
            DumpType,
            error::Result<Box<dyn std::io::Write + std::marker::Send>>,
        )> = vec![];

        match &self.pcon {
            None => (),
            Some(paths) => {
                for path in paths {
                    outputs.push((DumpType::Pcon, create(path)));
                }
            }
        }

        match &self.solid {
            None => (),
            Some(paths) => {
                for path in paths {
                    outputs.push((DumpType::Solid, create(path)));
                }
            }
        }

        match &self.csv {
            None => {
                if outputs.is_empty() {
                    outputs.push((
                        DumpType::Csv,
                        Ok(Box::new(std::io::BufWriter::new(std::io::stdout()))),
                    ))
                }
            }
            Some(paths) => {
                for path in paths {
                    outputs.push((DumpType::Csv, create(path)));
                }
            }
        }

        outputs
    }

    /// Get abundance
    pub fn abundance(&self) -> crate::CountTypeNoAtomic {
        self.abundance
    }
}

fn create<P>(path: P) -> error::Result<Box<dyn std::io::Write + std::marker::Send>>
where
    P: std::convert::AsRef<std::path::Path>,
{
    let file = std::fs::File::create(path)?;
    let buffer = std::io::BufWriter::new(file);
    let boxed = Box::new(buffer);

    Ok(boxed)
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::io::Write as _;

    #[cfg(not(feature = "parallel"))]
    #[test]
    fn basic() {
        let subcmd = Count {
            inputs: None,
            format: None,
            pcon: None,
            csv: None,
            solid: None,
            kmer_size: 32,
            abundance: Some(0),
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
            format: None,
            pcon: None,
            csv: None,
            solid: None,
            kmer_size: 32,
            abundance: None,
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

        let output = tempfile::NamedTempFile::new()?;

        let count = Count {
            inputs: Some(vec![
                input1.path().to_path_buf(),
                input2.path().to_path_buf(),
            ]),
            format: None,
            pcon: None,
            csv: None,
            solid: Some(vec![output.path().to_path_buf()]),
            kmer_size: 32,
            abundance: Some(2),
            record_buffer: Some(512),
        };

        let mut content = Vec::new();
        count.inputs()?.read_to_end(&mut content)?;
        assert_eq!(content, b">test\nTACG\n");

        assert_eq!(count.kmer_size(), 31);
        assert_eq!(count.abundance(), 2);
        assert_eq!(count.outputs()[0].0, DumpType::Solid);
        assert_eq!(count.record_buffer(), 512);

        let count = Count {
            inputs: Some(vec![
                input1.path().to_path_buf(),
                input2.path().to_path_buf(),
            ]),
            format: None,
            pcon: None,
            csv: None,
            solid: None,
            kmer_size: 32,
            abundance: Some(2),
            record_buffer: Some(512),
        };

        assert_eq!(count.outputs()[0].0, DumpType::Pcon);

        Ok(())
    }

    #[test]
    fn dump() -> error::Result<()> {
        let mut input1 = tempfile::NamedTempFile::new()?;
        input1.write_all(b">test\nATCG\n")?;

        let output = tempfile::NamedTempFile::new()?;

        let dump = Dump {
            input: Some(input1.path().to_path_buf()),
            pcon: None,
            csv: None,
            solid: Some(vec![output.path().to_path_buf()]),
            abundance: 2,
        };

        let mut content = Vec::new();
        dump.input()?.read_to_end(&mut content)?;
        assert_eq!(content, b">test\nATCG\n");

        assert_eq!(dump.abundance(), 2);
        assert_eq!(dump.outputs()[0].0, DumpType::Solid);

        Ok(())
    }
}
