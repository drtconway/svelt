use std::{iter::zip, sync::Arc};

use datafusion::{
    arrow::{
        array::{ArrayRef, ListBuilder, PrimitiveBuilder},
        datatypes::{DataType, Field, FieldRef, UInt64Type},
    },
    common::cast::{as_string_array, as_uint32_array},
    error::DataFusionError,
    logical_expr::{ColumnarValue, ScalarUDF, Volatility},
    prelude::SimpleScalarUDF,
};

use crate::{errors::SveltError, kmers::Kmer};

pub fn kmer_item() -> FieldRef {
    Arc::new(Field::new("kmer", DataType::UInt64, false))
}

pub fn kmer_list_type() -> DataType {
    DataType::List(kmer_item())
}

pub fn kmers_fwd() -> ScalarUDF {
    ScalarUDF::from(SimpleScalarUDF::new(
        "kmers_fwd",
        vec![DataType::UInt32, DataType::Utf8],
        kmer_list_type(),
        Volatility::Immutable,
        Arc::new(|args| do_kmers_fwd(args)),
    ))
}

fn do_kmers_fwd(args: &[ColumnarValue]) -> std::result::Result<ColumnarValue, DataFusionError> {
    let args = ColumnarValue::values_to_arrays(args)?;

    let ks = as_uint32_array(&args[0])?;
    let seqs = as_string_array(&args[1])?;

    let values_builder = PrimitiveBuilder::<UInt64Type>::new();
    let mut result_builder = ListBuilder::new(values_builder).with_field(kmer_item());
    for (k, seq) in zip(ks.iter(), seqs.iter()) {
        let k = if let Some(k) = k {
            k
        } else {
            return Err(DataFusionError::External(Box::new(
                SveltError::UnexpectedNull(String::from("kmers_fwd:k")),
            )));
        };
        let seq = if let Some(seq) = seq {
            seq
        } else {
            return Err(DataFusionError::External(Box::new(
                SveltError::UnexpectedNull(String::from("kmers_fwd:seq")),
            )));
        };
        let mut fwd = Vec::new();
        Kmer::with_many_both(k as usize, &seq, |x, _y| {
            fwd.push(Some(x.0));
        });
        result_builder.append_value(fwd);
    }

    let result_array = result_builder.finish();

    Ok(ColumnarValue::from(Arc::new(result_array) as ArrayRef))
}

pub fn kmers_rev() -> ScalarUDF {
    ScalarUDF::from(SimpleScalarUDF::new(
        "kmers_rev",
        vec![DataType::UInt32, DataType::Utf8],
        kmer_list_type(),
        Volatility::Immutable,
        Arc::new(|args| do_kmers_rev(args)),
    ))
}

fn do_kmers_rev(args: &[ColumnarValue]) -> std::result::Result<ColumnarValue, DataFusionError> {
    let args = ColumnarValue::values_to_arrays(args)?;

    let ks = as_uint32_array(&args[0])?;
    let seqs = as_string_array(&args[1])?;

    let values_builder = PrimitiveBuilder::<UInt64Type>::new();
    let mut result_builder = ListBuilder::new(values_builder).with_field(kmer_item());
    for (k, seq) in zip(ks.iter(), seqs.iter()) {
        let k = if let Some(k) = k {
            k
        } else {
            return Err(DataFusionError::External(Box::new(
                SveltError::UnexpectedNull(String::from("kmers_rev:k")),
            )));
        };
        let seq = if let Some(seq) = seq {
            seq
        } else {
            return Err(DataFusionError::External(Box::new(
                SveltError::UnexpectedNull(String::from("kmers_rev:seq")),
            )));
        };
        let mut rev = Vec::new();
        Kmer::with_many_both(k as usize, &seq, |_x, y| {
            rev.push(Some(y.0));
        });
        result_builder.append_value(rev);
    }

    let result_array = result_builder.finish();

    Ok(ColumnarValue::from(Arc::new(result_array) as ArrayRef))
}

#[cfg(test)]
mod tests {
    use datafusion::{
        arrow::{
            array::{Array, GenericListArray, RecordBatch, StringArray, UInt32Array, UInt64Array},
            datatypes::{Schema, SchemaRef},
        },
        prelude::{SessionContext, col},
    };

    use super::*;

    fn make_record(k: u32, seq: &str) -> RecordBatch {
        let schema = SchemaRef::new(Schema::new(vec![
            Field::new("k", DataType::UInt32, false),
            Field::new("seq", DataType::Utf8, true),
        ]));

        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt32Array::from(vec![k])),
                Arc::new(StringArray::from(vec![seq])),
            ],
        )
        .unwrap()
    }

    #[tokio::test]
    async fn test_kmers_fwd() {
        let f = kmers_fwd();
        let k: usize = 11;
        let seq = "GTCCCTCTGTCTCTGCCAACCAGTTAACCTGCTGCTTCCTGGAGGAA";
        let recs = make_record(k as u32, seq);
        let ctx = SessionContext::new();
        let df = ctx.read_batch(recs).unwrap();
        let df = df
            .with_column("kmers", f.call(vec![col("k"), col("seq")]))
            .unwrap();
        let res = df.collect().await.unwrap();
        assert_eq!(res.len(), 1);
        let cols = res[0].columns();
        let kmers = cols[2]
            .as_any()
            .downcast_ref::<GenericListArray<i32>>()
            .unwrap();
        assert_eq!(kmers.len(), 1);
        let xs = kmers.value(0);
        let xs = xs.as_any().downcast_ref::<UInt64Array>().unwrap();
        assert_eq!(xs.len(), seq.len() + 1 - k);
        assert_eq!(Kmer::from_u64(xs.value(0)).render(k), "GTCCCTCTGTC");
        assert_eq!(Kmer::from_u64(xs.value(11)).render(k), "TCTGCCAACCA");
        assert_eq!(Kmer::from_u64(xs.value(36)).render(k), "TCCTGGAGGAA");
    }

    #[tokio::test]
    async fn test_kmers_rev() {
        let f = kmers_rev();
        let k: usize = 11;
        let seq = "GTCCCTCTGTCTCTGCCAACCAGTTAACCTGCTGCTTCCTGGAGGAA";
        let recs = make_record(k as u32, seq);
        let ctx = SessionContext::new();
        let df = ctx.read_batch(recs).unwrap();
        let df = df
            .with_column("kmers", f.call(vec![col("k"), col("seq")]))
            .unwrap();
        let res = df.collect().await.unwrap();
        assert_eq!(res.len(), 1);
        let cols = res[0].columns();
        let kmers = cols[2]
            .as_any()
            .downcast_ref::<GenericListArray<i32>>()
            .unwrap();
        assert_eq!(kmers.len(), 1);
        let xs = kmers.value(0);
        let xs = xs.as_any().downcast_ref::<UInt64Array>().unwrap();
        assert_eq!(xs.len(), seq.len() + 1 - k);
        assert_eq!(Kmer::from_u64(xs.value(0)).render(k), "GACAGAGGGAC");
        assert_eq!(Kmer::from_u64(xs.value(11)).render(k), "TGGTTGGCAGA");
        assert_eq!(Kmer::from_u64(xs.value(36)).render(k), "TTCCTCCAGGA");
    }
}
