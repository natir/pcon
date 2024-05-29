//! Pcon utils function

/* std use */

/* crate use */

/* project use */

/// Initialize counter
pub fn init_data<T>(k: u8, value: T) -> Box<[T]>
where
    T: std::marker::Sized + std::clone::Clone,
{
    vec![value; cocktail::kmer::get_hash_space_size(k) as usize].into_boxed_slice()
}

#[cfg(feature = "parallel")]
/// Perform transmutation on box
pub fn transmute<I, O>(data: &[I]) -> &[O]
where
    I: std::marker::Sized,
    O: std::marker::Sized,
{
    unsafe { std::mem::transmute::<&[I], &[O]>(data) }
}

#[cfg(feature = "parallel")]
/// Perform transmutation on box
pub fn transmute_box<I, O>(data: Box<[I]>) -> Box<[O]>
where
    I: std::marker::Sized,
    O: std::marker::Sized,
{
    unsafe { std::mem::transmute::<Box<[I]>, Box<[O]>>(data) }
}

#[cfg(feature = "parallel")]
/// Populate record buffer with content of iterator
pub fn populate_buffer(
    iter: &mut noodles::fasta::reader::Records<'_, Box<dyn std::io::BufRead>>,
    records: &mut Vec<noodles::fasta::Record>,
    record_buffer: u64,
) -> bool {
    records.clear();

    for i in 0..record_buffer {
        if let Some(Ok(record)) = iter.next() {
            records.push(record);
        } else {
            records.truncate(i as usize);
            return false;
        }
    }

    true
}

#[cfg(all(feature = "parallel", feature = "fastq"))]
/// Populate record buffer with content of iterator
pub fn populate_bufferq(
    iter: &mut noodles::fastq::reader::Records<'_, Box<dyn std::io::BufRead>>,
    records: &mut Vec<noodles::fastq::Record>,
    record_buffer: u64,
) -> bool {
    records.clear();

    for i in 0..record_buffer {
        if let Some(Ok(record)) = iter.next() {
            records.push(record);
        } else {
            records.truncate(i as usize);
            return false;
        }
    }

    true
}
