use std::io::{BufReader, Error, ErrorKind};

use autocompress::autodetect_open;

use super::{InnerReader, SequenceIterator};

pub struct FastaSequenceIterator {
    pub(crate) reader: noodles::fasta::io::Reader<InnerReader>,
}

impl FastaSequenceIterator {
    pub fn new(filename: &str) -> std::io::Result<FastaSequenceIterator> {
        let reader = autodetect_open(filename)?;
        let reader = BufReader::new(reader);
        let reader = noodles::fasta::io::reader::Builder::default().build_from_reader(reader)?;
        Ok(FastaSequenceIterator { reader })
    }

    pub(crate) fn read_one(&mut self) -> std::io::Result<Option<(String, String)>> {
        let mut definition = String::new();
        let mut sequence = Vec::new();

        let r1 = self.reader.read_definition(&mut definition)?;
        let definition = String::from(&definition[1..]);

        if r1 == 0 {
            return Ok(None);
        }

        self.reader.read_sequence(&mut sequence)?;
        let sequence = String::from_utf8(sequence).map_err(|e| Error::new(ErrorKind::Other, e))?;

        Ok(Some((definition, sequence)))
    }
}

impl Iterator for FastaSequenceIterator {
    type Item = std::io::Result<(String, String)>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.read_one() {
            Ok(Some(item)) => Some(Ok(item)),
            Ok(None) => None,
            Err(error) => Some(Err(error)),
        }
    }
}

impl SequenceIterator for FastaSequenceIterator {}
