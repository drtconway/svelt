use std::{collections::HashMap, rc::Rc};

use noodles::vcf::{Header, Record};

use crate::{chroms::ChromSet, vcf_reader::VcfReader};



pub struct RecordSeeker {
    reader: VcfReader,
    header: Rc<Header>,
    number_read: u32,
    buffer: HashMap<u32, Record>  
}

impl RecordSeeker {
    pub fn new(path: &str, chroms: Rc<ChromSet>) -> std::io::Result<RecordSeeker> {

        let reader = VcfReader::new(path, chroms)?;
        let header = Rc::new(reader.header.clone());
        let buffer = HashMap::new();
        Ok(RecordSeeker { reader, header, number_read: 0, buffer })
    }

    pub fn take(&mut self, rn: u32) -> std::io::Result<Option<(Rc<Header>, Record)>> {
        while self.number_read <= rn {
            let mut rec = Record::default();
            let res = self.reader.reader.read_record(&mut rec)?;
            if res == 0 {
                // EOF
                break;
            }
            self.buffer.insert(self.number_read, rec);
            self.number_read += 1;
        }
        Ok(self.buffer.remove(&rn).map(|rec| (self.header.clone(), rec)))
    }
}