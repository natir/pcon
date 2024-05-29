/* std use */
use std::io::Read as _;

/* 3rd party use */
use biotest::Format as _;

/* local use */
pub mod constant;

mod count {
    /* local use */
    use super::*;

    #[cfg(not(any(feature = "count_u16", feature = "count_u32", feature = "count_u64")))]
    #[test]
    fn from_stdin_to_stdout() -> anyhow::Result<()> {
        let mut rng = biotest::rand();
        let generator = biotest::Fasta::builder().sequence_len(150).build()?;

        let mut buffer = Vec::new();
        generator.records(&mut buffer, &mut rng, 100)?;

        let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
        cmd.args(["count", "-k", "5"]).write_stdin(buffer);

        let assert = cmd.assert();

        assert
            .success()
            .stderr(b"" as &[u8])
            .stdout(constant::TRUTH_PCON);
        Ok(())
    }

    #[cfg(not(any(feature = "count_u16", feature = "count_u32", feature = "count_u64")))]
    #[test]
    fn from_file_to_stdout() -> anyhow::Result<()> {
        let mut rng = biotest::rand();
        let generator = biotest::Fasta::builder().sequence_len(150).build()?;

        let input_temp = tempfile::NamedTempFile::new()?;
        let input_path = input_temp.path();

        generator.create(input_path, &mut rng, 100)?;

        let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
        cmd.args([
            "count",
            "-k",
            "5",
            "-i",
            &format!("{}", input_path.display()),
        ]);

        let assert = cmd.assert();

        assert
            .success()
            .stderr(b"" as &[u8])
            .stdout(constant::TRUTH_PCON);

        Ok(())
    }

    #[cfg(not(any(feature = "count_u16", feature = "count_u32", feature = "count_u64")))]
    #[test]
    fn from_file_to_file() -> anyhow::Result<()> {
        let mut rng = biotest::rand();
        let generator = biotest::Fasta::builder().sequence_len(150).build()?;

        let input_temp = tempfile::NamedTempFile::new()?;
        let input_path = input_temp.path();

        generator.create(input_path, &mut rng, 100)?;

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

    #[cfg(not(any(feature = "count_u16", feature = "count_u32", feature = "count_u64")))]
    #[test]
    fn count_to_csv() -> anyhow::Result<()> {
        let mut rng = biotest::rand();
        let generator = biotest::Fasta::builder().sequence_len(150).build()?;

        let mut buffer = Vec::new();
        generator.records(&mut buffer, &mut rng, 100)?;

        let mut output_temp = tempfile::NamedTempFile::new()?;
        let output_path = output_temp.path();

        let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
        cmd.args([
            "count",
            "-k",
            "5",
            "-a",
            "35",
            "-c",
            &format!("{}", output_path.display()),
        ])
        .write_stdin(buffer);

        let assert = cmd.assert();

        assert.success().stderr(b"" as &[u8]).stdout(b"" as &[u8]);

        let mut output = vec![];
        output_temp.read_to_end(&mut output)?;
        assert_eq!(output, constant::TRUTH_CSV);

        Ok(())
    }

    #[cfg(not(any(feature = "count_u16", feature = "count_u32", feature = "count_u64")))]
    #[test]
    fn count_to_solid() -> anyhow::Result<()> {
        let mut rng = biotest::rand();
        let generator = biotest::Fasta::builder().sequence_len(150).build()?;

        let mut buffer = Vec::new();
        generator.records(&mut buffer, &mut rng, 100)?;

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
        .write_stdin(buffer);

        let assert = cmd.assert();

        println!("{:?}", assert.get_output().stdout);

        assert.success().stderr(b"" as &[u8]).stdout(b"" as &[u8]);

        let mut output = vec![];
        output_temp.read_to_end(&mut output)?;
        assert_eq!(output, constant::TRUTH_SOLID);

        Ok(())
    }
}
