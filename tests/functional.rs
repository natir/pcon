/* std use */
use std::io::Read as _;
use std::io::Write as _;

/* 3rd party use */
use rand::prelude::*;

/* integrate constant */
pub mod constant;

/* Utils function */
fn fasta_buffer(nb_reads: usize, seq_len: usize) -> Vec<u8> {
    let mut ret = Vec::new();

    let mut rng = rand::rngs::StdRng::from_seed(constant::SEED);
    let nucs = [b'A', b'C', b'T', b'G', b'a', b'c', b't', b'g'];

    for _ in 0..nb_reads {
        let seq = (0..seq_len)
            .map(|_| *nucs.choose(&mut rng).unwrap())
            .collect::<Vec<u8>>();

        ret.extend(b">random\n");
        ret.extend(&seq[..]);
        ret.extend(b"\n");
    }

    ret
}

fn write_fasta<P>(path: P) -> std::io::Result<()>
where
    P: std::convert::AsRef<std::path::Path>,
{
    let mut output = std::fs::File::create(&path)?;

    output.write_all(&fasta_buffer(100, 150))
}

#[cfg(test)]
mod tests {

    use super::*;

    mod count {
        use super::*;

        #[test]
        fn from_stdin_to_stdout() -> std::io::Result<()> {
            let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
            cmd.args(["count", "-k", "5"])
                .write_stdin(fasta_buffer(100, 150));

            let assert = cmd.assert();

            assert
                .success()
                .stderr(b"" as &[u8])
                .stdout(constant::TRUTH_PCON);
            Ok(())
        }

        #[test]
        fn from_file_to_stdout() -> std::io::Result<()> {
            let input_temp = tempfile::NamedTempFile::new()?;
            let input_path = input_temp.path();

            write_fasta(input_path)?;

            let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
            cmd.args([
                "count",
                "-k",
                "5",
                "-i",
                &format!("{}", input_path.display()),
            ]);

            let assert = cmd.assert();

            println!("{:?}", assert.get_output().stdout);

            assert
                .success()
                .stderr(b"" as &[u8])
                .stdout(constant::TRUTH_PCON);

            Ok(())
        }

        #[test]
        fn from_file_to_file() -> std::io::Result<()> {
            let input_temp = tempfile::NamedTempFile::new()?;
            let input_path = input_temp.path();

            write_fasta(input_path)?;

            let mut output_temp = tempfile::NamedTempFile::new()?;
            let output_path = output_temp.path();

            let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
            cmd.args([
                "count",
                "-k",
                "5",
                "-i",
                &format!("{}", input_path.display()),
                "-p",
                &format!("{}", output_path.display()),
            ]);

            let assert = cmd.assert();

            assert.success().stderr(b"" as &[u8]).stdout(b"" as &[u8]);

            let mut output = vec![];
            output_temp.read_to_end(&mut output)?;
            assert_eq!(output, constant::TRUTH_PCON);

            Ok(())
        }
    }

    mod dump {
        use super::*;

        #[test]
        fn from_stdin_to_stdout() -> std::io::Result<()> {
            let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
            cmd.args(["dump", "-a", "1"])
                .write_stdin(constant::TRUTH_PCON);

            let assert = cmd.assert();

            assert
                .success()
                .stderr(b"" as &[u8])
                .stdout(constant::TRUTH_CSV);
            Ok(())
        }

        #[test]
        fn from_file_to_stdout() -> std::io::Result<()> {
            let mut input_temp = tempfile::NamedTempFile::new()?;
            {
                input_temp.write_all(constant::TRUTH_PCON)?;
            }
            let input_path = input_temp.path();

            let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
            cmd.args([
                "dump",
                "-a",
                "1",
                "-i",
                &format!("{}", input_path.display()),
            ]);

            let assert = cmd.assert();

            assert
                .success()
                .stderr(b"" as &[u8])
                .stdout(constant::TRUTH_CSV);

            Ok(())
        }

        #[test]
        fn from_file_to_file() -> std::io::Result<()> {
            let mut input_temp = tempfile::NamedTempFile::new()?;
            {
                input_temp.write_all(constant::TRUTH_PCON)?;
            }
            let input_path = input_temp.path();

            let mut output_temp = tempfile::NamedTempFile::new()?;
            let output_path = output_temp.path();

            let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
            cmd.args([
                "dump",
                "-a",
                "1",
                "-i",
                &format!("{}", input_path.display()),
                "-c",
                &format!("{}", output_path.display()),
            ]);

            let assert = cmd.assert();

            assert.success().stderr(b"" as &[u8]).stdout(b"" as &[u8]);

            let mut output = vec![];
            output_temp.read_to_end(&mut output)?;
            assert_eq!(output, constant::TRUTH_CSV);

            Ok(())
        }
    }

    mod other_output {
        use super::*;

        #[test]
        fn count_csv() -> std::io::Result<()> {
            let mut output_temp = tempfile::NamedTempFile::new()?;
            let output_path = output_temp.path();

            let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
            cmd.args([
                "count",
                "-k",
                "5",
                "-c",
                &format!("{}", output_path.display()),
            ])
            .write_stdin(fasta_buffer(100, 150));

            let assert = cmd.assert();

            assert.success().stderr(b"" as &[u8]).stdout(b"" as &[u8]);

            let mut output = vec![];
            output_temp.read_to_end(&mut output)?;
            assert_eq!(output, constant::TRUTH_CSV);

            Ok(())
        }

        #[test]
        fn count_solid() -> std::io::Result<()> {
            let mut output_temp = tempfile::NamedTempFile::new()?;
            let output_path = output_temp.path();

            let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
            cmd.args([
                "count",
                "-k",
                "5",
                "-s",
                &format!("{}", output_path.display()),
            ])
            .write_stdin(fasta_buffer(100, 150));

            let assert = cmd.assert();

            println!("{:?}", assert.get_output().stdout);

            assert.success().stderr(b"" as &[u8]).stdout(b"" as &[u8]);

            let mut output = vec![];
            output_temp.read_to_end(&mut output)?;
            assert_eq!(output, constant::TRUTH_SOLID);

            Ok(())
        }

        #[test]
        fn dump_pcon() -> std::io::Result<()> {
            let mut output_temp = tempfile::NamedTempFile::new()?;
            let output_path = output_temp.path();

            let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
            cmd.args([
                "dump",
                "-a",
                "1",
                "-p",
                &format!("{}", output_path.display()),
            ])
            .write_stdin(constant::TRUTH_PCON);

            let assert = cmd.assert();

            assert.success().stderr(b"" as &[u8]).stdout(b"" as &[u8]);

            let mut output = vec![];
            output_temp.read_to_end(&mut output)?;
            assert_eq!(output, constant::TRUTH_PCON);

            Ok(())
        }

        #[test]
        fn dump_solid() -> std::io::Result<()> {
            let mut output_temp = tempfile::NamedTempFile::new()?;
            let output_path = output_temp.path();

            let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
            cmd.args([
                "dump",
                "-a",
                "1",
                "-s",
                &format!("{}", output_path.display()),
            ])
            .write_stdin(constant::TRUTH_PCON);

            let assert = cmd.assert();

            assert.success().stderr(b"" as &[u8]).stdout(b"" as &[u8]);

            let mut output = vec![];
            output_temp.read_to_end(&mut output)?;
            assert_eq!(output, constant::TRUTH_SOLID);

            Ok(())
        }
    }
}
