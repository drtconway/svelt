use std::io::{BufReader, };

use autocompress::{Processor, io::ProcessorReader};


pub trait SequenceIterator: Iterator<Item = std::io::Result<(String, String)>> {}

type InnerReader = BufReader<
    ProcessorReader<Box<dyn Processor + Send + Unpin + 'static>, BufReader<std::fs::File>>,
>;

pub mod fasta;

pub mod vcf;
