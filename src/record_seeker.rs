use std::{collections::HashMap, rc::Rc};

use noodles::vcf::{Header, Record};

use crate::{chroms::ChromSet, vcf_reader::VcfReader};



pub struct RecordSeeker {
    header: Rc<Header>,
    buffer: HashMap<u32, Record>  
}

impl RecordSeeker {
    pub fn new(path: &str, chroms: Rc<ChromSet>) -> std::io::Result<RecordSeeker> {
        let mut buffer = HashMap::new();

        let mut reader = VcfReader::new(path, chroms)?;
        let header = Rc::new(reader.header.clone());

        for (rn, rec) in reader.reader.records().enumerate() {
            let rec = rec?;
            buffer.insert(rn as u32, rec);
        }
        log::info!("read {} records from '{}'", buffer.len(), path);

        Ok(RecordSeeker { header, buffer })
    }

    pub fn take(&mut self, rn: u32) -> Option<(Rc<Header>, Record)> {
        self.buffer.remove(&rn).map(|rec| (self.header.clone(), rec))
    }
}