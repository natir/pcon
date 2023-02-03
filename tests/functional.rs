/* std use */
use std::io::Read as _;
use std::io::Seek as _;
use std::io::Write as _;

/* 3rd party use */
use rand::prelude::*;

fn write_fasta<P>(path: P) -> std::io::Result<()>
where
    P: std::convert::AsRef<std::path::Path>,
{
    let mut output = std::fs::File::create(&path)?;

    let mut rng = rand::rngs::StdRng::from_seed([42; 32]);
    let nucs = [b'A', b'C', b'T', b'G', b'a', b'c', b't', b'g'];

    for _ in 0..100 {
        let seq = (0..150)
            .map(|_| *nucs.choose(&mut rng).unwrap())
            .collect::<Vec<u8>>();

        output.write_all(b">random")?;
        output.write_all(&seq[..])?;
    }

    Ok(())
}

fn run_finish(mut child: std::process::Child) -> std::io::Result<Vec<u8>> {
    if !child.wait()?.success() {
        let mut stdout = String::new();
        let mut stderr = String::new();

        child
            .stdout
            .take()
            .expect("No stdout")
            .read_to_string(&mut stdout)?;
        child
            .stderr
            .take()
            .expect("No stderr")
            .read_to_string(&mut stderr)?;

        println!("stdout: {stdout}");
        println!("stderr: {stderr}");

        Err(std::io::Error::new(
            std::io::ErrorKind::Other,
            "Command not finish",
        ))
    } else {
        let mut stdout_content = Vec::new();
        let mut stdout = child.stdout.take().expect("No stdout");
        stdout.read_to_end(&mut stdout_content)?;

        Ok(stdout_content)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn count_to_stdout() -> std::io::Result<()> {
        let tempfile = tempfile::NamedTempFile::new()?;
        let path = tempfile.path();

        write_fasta(path)?;

        let child = std::process::Command::new("./target/debug/pcon")
            .args(&["count", "-k", "5", "-i", &format!("{}", path.display())])
            .stderr(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .spawn()?;

        let stdout = run_finish(child)?;

        assert_eq!(
            stdout,
            vec![
                5, 1, 31, 139, 8, 0, 0, 0, 0, 0, 4, 255, 237, 208, 129, 0, 0, 0, 0, 128, 160, 253,
                169, 7, 249, 16, 98, 233, 1, 1, 120, 117, 170, 178, 0, 2, 0, 0
            ]
        );

        Ok(())
    }

    #[test]
    fn count_from_stdin() -> std::io::Result<()> {
        let tempfile = tempfile::NamedTempFile::new()?;
        let path = tempfile.path();

        write_fasta(path)?;

        let child = std::process::Command::new("./target/debug/pcon")
            .args(&["count", "-k", "5"])
            .stdin(std::fs::File::open(path)?)
            .stderr(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .spawn()?;

        let stdout = run_finish(child)?;

        assert_eq!(
            stdout,
            vec![
                5, 1, 31, 139, 8, 0, 0, 0, 0, 0, 4, 255, 237, 208, 129, 0, 0, 0, 0, 128, 160, 253,
                169, 7, 249, 16, 98, 233, 1, 1, 120, 117, 170, 178, 0, 2, 0, 0
            ]
        );

        Ok(())
    }
}
