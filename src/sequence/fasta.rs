use std::io::{BufReader, Error, ErrorKind};

use noodles::fasta::record::Definition;

use super::{InnerReader, SequenceIterator};

pub struct FastaSequenceIterator {
    pub(crate) reader: noodles::fasta::io::Reader<InnerReader>,
}

impl FastaSequenceIterator {
    pub fn new(filename: &str) -> std::io::Result<FastaSequenceIterator> {
        let (reader, _) = niffler::send::from_path(std::path::Path::new(filename))
            .map_err(|e| Error::new(ErrorKind::Other, e))?;
        let reader = BufReader::new(reader);
        let reader = noodles::fasta::io::reader::Builder::default().build_from_reader(reader)?;
        Ok(FastaSequenceIterator { reader })
    }

    pub(crate) fn read_one(&mut self) -> std::io::Result<Option<(String, String)>> {
        let mut definition = Definition::default();
        let mut sequence = Vec::new();

        let r1 = self.reader.read_definition(&mut definition)?;

        if r1 == 0 {
            return Ok(None);
        }

        let label = definition.name().to_string();

        self.reader.read_sequence(&mut sequence)?;
        let sequence = String::from_utf8(sequence).map_err(|e| Error::new(ErrorKind::Other, e))?;

        Ok(Some((label, sequence)))
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
