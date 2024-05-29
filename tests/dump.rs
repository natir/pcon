/* std use */
use std::io::Read as _;
use std::io::Write as _;

/* 3rd party use */

/* local use */
pub mod constant;

mod dump {
    /* local use */
    use super::*;

    #[cfg(not(any(feature = "count_u16", feature = "count_u32", feature = "count_u64")))]
    #[test]
    fn from_stdin_to_stdout() -> std::io::Result<()> {
        let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
        cmd.args(["dump", "-a", "35"])
            .write_stdin(constant::TRUTH_PCON);

        let assert = cmd.assert();

        assert
            .success()
            .stderr(b"" as &[u8])
            .stdout(constant::TRUTH_CSV);
        Ok(())
    }

    #[cfg(not(any(feature = "count_u16", feature = "count_u32", feature = "count_u64")))]
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
            "35",
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

    #[cfg(not(any(feature = "count_u16", feature = "count_u32", feature = "count_u64")))]
    #[test]
    fn from_file_to_file() -> std::io::Result<()> {
        let mut input_temp = tempfile::NamedTempFile::new()?;
        input_temp.write_all(constant::TRUTH_PCON)?;
        let input_path = input_temp.path();

        let mut output_temp = tempfile::NamedTempFile::new()?;
        let output_path = output_temp.path();

        let mut cmd = assert_cmd::Command::cargo_bin("pcon").unwrap();
        cmd.args([
            "dump",
            "-a",
            "35",
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

    #[cfg(not(any(feature = "count_u16", feature = "count_u32", feature = "count_u64")))]
    #[test]
    fn dump_to_pcon() -> anyhow::Result<()> {
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

    #[cfg(not(any(feature = "count_u16", feature = "count_u32", feature = "count_u64")))]
    #[test]
    fn dump_to_solid() -> anyhow::Result<()> {
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
