/* std use */
use std::io::Read as _;

/* 3rd party use */
use biotest::Format as _;

/* local use */
pub mod constant;

mod minicount {
    /* local use */
    use super::*;

    #[test]
    fn from_stdin_to_stdout() -> anyhow::Result<()> {
        let mut rng = biotest::rand();
        let generator = biotest::Fasta::builder()
            .sequence_len(150)
            .sequence(biotest::values::Nucleotides::DnaUpper)
            .build()?;

        let mut buffer = Vec::new();
        generator.records(&mut buffer, &mut rng, 100)?;

        let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
        cmd.args(["mini-count", "-k", "5", "-m", "3", "-a", "20"])
            .write_stdin(buffer);

        let assert = cmd.assert();

        let mut output = assert.get_output().stdout.to_vec();
        output.sort_unstable();

        let mut truth = constant::TRUTH_CSV_MINI.to_vec();
        truth.sort_unstable();

        assert.success().stderr(b"" as &[u8]);

        Ok(())
    }

    #[cfg(not(feature = "parallel"))]
    #[ignore]
    #[test]
    fn from_file_to_stdout() -> anyhow::Result<()> {
        let mut rng = biotest::rand();
        let generator = biotest::Fasta::builder().sequence_len(150).build()?;

        let input_temp = tempfile::NamedTempFile::new()?;
        let input_path = input_temp.path();

        generator.create(input_path, &mut rng, 100)?;

        let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
        cmd.args([
            "mini-count",
            "-k",
            "5",
            "-m",
            "3",
            "-a",
            "20",
            "-i",
            &format!("{}", input_path.display()),
        ]);

        let assert = cmd.assert();

        let mut output = assert.get_output().stdout.to_vec();
        output.sort_unstable();

        let mut truth = constant::TRUTH_CSV_MINI.to_vec();
        truth.sort_unstable();

        assert_eq!(output, truth);

        Ok(())
    }

    #[test]
    fn from_file_to_file() -> anyhow::Result<()> {
        let mut rng = biotest::rand();
        let generator = biotest::Fasta::builder().sequence_len(150).build()?;

        let input_temp = tempfile::NamedTempFile::new()?;
        let input_path = input_temp.path();

        generator.create(input_path, &mut rng, 100)?;

        let output_temp = tempfile::NamedTempFile::new()?;
        let output_path = output_temp.path();

        let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
        cmd.args([
            "mini-count",
            "-k",
            "5",
            "-m",
            "3",
            "-a",
            "20",
            "-i",
            &format!("{}", input_path.display()),
            "-c",
            &format!("{}", output_path.display()),
        ]);

        let assert = cmd.assert();

        let mut output = Vec::new();
        std::fs::File::open(output_path)?.read_to_end(&mut output)?;
        output.sort_unstable();

        let mut truth = constant::TRUTH_CSV_MINI.to_vec();
        truth.sort_unstable();

        assert.success().stderr(b"" as &[u8]);

        Ok(())
    }
}
