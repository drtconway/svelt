use noodles::vcf::{Header, variant::Record};

use crate::errors::{Context, VariantContext, as_io_error};

/// Traverse a VCF record to force any errors to become apparent.
pub fn traverse_record<R: Record>(rec: &R, header: &Header) -> std::io::Result<()> {
    let chrom = rec.reference_sequence_name(header)?;

    let position = match rec.variant_start() {
        Some(variant_start) => variant_start?.get(),
        None => 0,
    };

    VariantContext::new(chrom, position).with(|| {
        for id in rec.ids().iter() {
            let _ = id;
        }

        for base in rec.reference_bases().iter() {
            let _ = base?;
        }

        for alt in rec.alternate_bases().iter() {
            let _ = alt?;
        }

        match rec.quality_score() {
            Some(quality_score) => {
                let _ = quality_score?;
            }
            None => {}
        }

        for filter in rec.filters().iter(header) {
            let _ = filter?;
        }

        for info in rec.info().iter(header) {
            let info = info?;
            InfoContext::new(info.0).with(|| {
                match info.1 {
                    Some(value) => {
                        if let noodles::vcf::variant::record::info::field::Value::Array(array) = value {
                            match array {
                                noodles::vcf::variant::record::info::field::value::Array::Integer(
                                    values,
                                ) => {
                                    for value in values.iter() {
                                        let _ = value?;
                                    }
                                }
                                noodles::vcf::variant::record::info::field::value::Array::Float(
                                    values,
                                ) => {
                                    for value in values.iter() {
                                        let _ = value?;
                                    }
                                }
                                noodles::vcf::variant::record::info::field::value::Array::Character(
                                    values,
                                ) => {
                                    for value in values.iter() {
                                        let _ = value?;
                                    }
                                }
                                noodles::vcf::variant::record::info::field::value::Array::String(
                                    values,
                                ) => {
                                    for value in values.iter() {
                                        let _ = value?;
                                    }
                                }
                            }
                        }
                    }
                    None => {}
                }
                Ok(())
            })?;
        }

        let sample_names = header.sample_names();
        let samples = rec.samples()?;
        for (i, sample) in samples.iter().enumerate() {
            let sample_name = sample_names.get_index(i).unwrap();
            SampleContext::new(&sample_name).with(|| {
                for geno in sample.iter(header) {
                    let geno = geno?;
                    FmtContext::new(geno.0).with(|| {
                        match geno.1 {
                            Some(value) => {
                                match value {
                                    noodles::vcf::variant::record::samples::series::Value::Genotype(genotype) => {
                                        for allele in genotype.iter() {
                                            let _ = allele?;
                                        }
                                    },
                                    noodles::vcf::variant::record::samples::series::Value::Array(array) => {
                                        match array {
                                            noodles::vcf::variant::record::samples::series::value::Array::Integer(values) => {
                                                for value in values.iter() {
                                                    let _ = value?;
                                                }
                                            },
                                            noodles::vcf::variant::record::samples::series::value::Array::Float(values) => {
                                                for value in values.iter() {
                                                    let _ = value?;
                                                }
                                            },
                                            noodles::vcf::variant::record::samples::series::value::Array::Character(values) =>{
                                                for value in values.iter() {
                                                    let _ = value?;
                                                }
                                            },
                                            noodles::vcf::variant::record::samples::series::value::Array::String(values) => {
                                                for value in values.iter() {
                                                    let _ = value?;
                                                }
                                            },
                                        }
                                    },
                                    _ => {}
                                }
                            },
                            None => {},
                        }
                        Ok(())
                    })?;
                }
                Ok(())
            })?;
        }
        Ok(())
    })
}

pub struct InfoContext<'a> {
    name: &'a str,
}

impl<'a> InfoContext<'a> {
    pub fn new(name: &'a str) -> Self {
        InfoContext { name }
    }
}

impl<'a> Context for InfoContext<'a> {
    fn with<R, F: FnOnce() -> std::io::Result<R>>(&self, inner: F) -> std::io::Result<R> {
        (inner)().map_err(|e| {
            as_io_error(crate::errors::SveltError::BadInfoField(
                self.name.to_string(),
                Box::new(e),
            ))
        })
    }
}

pub struct SampleContext<'a> {
    name: &'a str,
}

impl<'a> SampleContext<'a> {
    pub fn new(name: &'a str) -> Self {
        SampleContext { name }
    }
}

impl<'a> Context for SampleContext<'a> {
    fn with<R, F: FnOnce() -> std::io::Result<R>>(&self, inner: F) -> std::io::Result<R> {
        (inner)().map_err(|e| {
            as_io_error(crate::errors::SveltError::BadSample(
                self.name.to_string(),
                Box::new(e),
            ))
        })
    }
}

pub struct FmtContext<'a> {
    name: &'a str,
}

impl<'a> FmtContext<'a> {
    pub fn new(name: &'a str) -> Self {
        FmtContext { name }
    }
}

impl<'a> Context for FmtContext<'a> {
    fn with<R, F: FnOnce() -> std::io::Result<R>>(&self, inner: F) -> std::io::Result<R> {
        (inner)().map_err(|e| {
            as_io_error(crate::errors::SveltError::BadFormatField(
                self.name.to_string(),
                Box::new(e),
            ))
        })
    }
}
