use std::{io::{Error, ErrorKind}, sync::Arc};

use datafusion::{
    arrow::{
        array::{
            Array, DictionaryArray, GenericByteArray, Int64Array, PrimitiveArray, PrimitiveBuilder,
            RecordBatch, StringArray, UInt16Array, UInt32Array,
        },
        datatypes::{DataType, Field, GenericStringType, Int64Type, Schema, UInt16Type, UInt32Type, UInt8Type},
    },
    prelude::{col, lit, DataFrame},
};

use crate::{
    heap::{Heap, HeapItem},
    row_key::RowKey,
};

pub async fn near_merge_table(tbl: DataFrame) -> std::io::Result<RecordBatch> {
    let tbl = tbl.filter(col("kind").not_eq(lit("BND")))?.sort_by(vec![
        col("kind"),
        col("chrom_id"),
        col("start"),
        col("end"),
    ])?;
    let batch = tbl.collect().await?;

    let mut lhs_row_id_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut rhs_row_id_builder = PrimitiveBuilder::<UInt32Type>::new();

    let mut lhs_heap: Heap<Row<'_>> = Heap::new();
    let mut lhs_itr = batch.iter().flat_map(|recs| MergeIterator::new(recs));
    let mut lhs = lhs_itr.next();

    let mut rhs_heap: Heap<Row<'_>> = Heap::new();
    let mut rhs_itr = batch.iter().flat_map(|recs| MergeIterator::new(recs));
    let mut rhs = rhs_itr.next();

    while let (Some(lhs_val), Some(rhs_val)) = (&lhs, &rhs) {
        //log::info!("iterating at {:?} and {:?}", lhs_val, rhs_val);

        // Kind
        if lhs_val.kind < rhs_val.kind {
            //log::info!("skip lhs: {}", lhs_val.kind);
            lhs_heap.clear();
            rhs_heap.clear();
            lhs = lhs_itr.next();
            continue;
        }
        if rhs_val.kind < lhs_val.kind {
            //log::info!("skip rhs: {}", rhs_val.kind);
            lhs_heap.clear();
            rhs_heap.clear();
            rhs = rhs_itr.next();
            continue;
        }

        // Chrom
        if lhs_val.chrom_id < rhs_val.chrom_id {
            //log::info!("skip lhs: {}", lhs_val.chrom_id);
            lhs_heap.clear();
            rhs_heap.clear();
            lhs = lhs_itr.next();
            continue;
        }
        if rhs_val.chrom_id < lhs_val.chrom_id {
            //log::info!("skip rhs: {}", rhs_val.chrom_id);
            lhs_heap.clear();
            rhs_heap.clear();
            rhs = rhs_itr.next();
            continue;
        }

        let lhs_row_key = RowKey::decode(lhs_val.row_id as u32);
        let lhs_start = lhs_val.start;
        let rhs_row_key = RowKey::decode(rhs_val.row_id as u32);
        let rhs_start = rhs_val.start;

        if lhs_start <= rhs_start {
            while let Some(rhs_item) = rhs_heap.front() {
                if rhs_item.start + 25 < lhs_start {
                    rhs_heap.pop();
                } else {
                    break;
                }
            }

            for rhs_item in rhs_heap.iter() {
                assert_eq!(lhs_val.kind, rhs_item.kind);
                assert_eq!(lhs_val.chrom_id, rhs_item.chrom_id);
                assert!(rhs_item.start + 25 >= lhs_start);
                if ((lhs_val.end as i64) - (rhs_item.end as i64)).abs() > 25 {
                    continue;
                }
                if lhs_val.row_id >= rhs_item.row_id {
                    continue;
                }
                let rhs_item_row_key = RowKey::decode(rhs_item.row_id as u32);
                if lhs_row_key.0 == rhs_item_row_key.0 {
                    continue;
                }
                lhs_row_id_builder.append_value(lhs_val.row_id as u32);
                rhs_row_id_builder.append_value(rhs_item.row_id as u32);
                //log::info!("joining-lhs {:?} and {:?}", lhs_val, rhs_item);
            }

            lhs_heap.push(lhs.take().unwrap());
            lhs = lhs_itr.next();
        } else {
            while let Some(lhs_item) = lhs_heap.front() {
                if lhs_item.start + 25 < rhs_start {
                    lhs_heap.pop();
                } else {
                    break;
                }
            }

            for lhs_item in lhs_heap.iter() {
                assert_eq!(rhs_val.kind, lhs_item.kind);
                assert_eq!(rhs_val.chrom_id, lhs_item.chrom_id);
                assert!(lhs_item.start + 25 >= rhs_start);
                if ((rhs_val.end as i64) - (lhs_item.end as i64)).abs() > 25 {
                    continue;
                }
                if rhs_val.row_id <= lhs_item.row_id {
                    continue;
                }
                let lhs_item_row_key = RowKey::decode(lhs_item.row_id as u32);
                if rhs_row_key.0 == lhs_item_row_key.0 {
                    continue;
                }
                lhs_row_id_builder.append_value(lhs_item.row_id as u32);
                rhs_row_id_builder.append_value(rhs_val.row_id as u32);
                //log::info!("joining-rhs {:?} and {:?}", rhs_val, lhs_item);
            }

            rhs_heap.push(rhs.take().unwrap());
            rhs = rhs_itr.next();
        }
    }

    let lhs_row_id_array = lhs_row_id_builder.finish();
    let rhs_row_id_array = rhs_row_id_builder.finish();

    let schema = Arc::new(Schema::new(vec![
        Field::new("lhs_row_id", DataType::UInt32, false),
        Field::new("rhs_row_id", DataType::UInt32, false),
    ]));

    let recs = RecordBatch::try_new(
        schema,
        vec![Arc::new(lhs_row_id_array), Arc::new(rhs_row_id_array)],
    )
    .map_err(|e| Error::new(ErrorKind::Other, e))?;

    Ok(recs)
}

struct MergeIterator<'a> {
    kind: &'a DictionaryArray<UInt8Type>,
    kind_values: &'a GenericByteArray<GenericStringType<i32>>,
    chrom_id: &'a PrimitiveArray<UInt16Type>,
    start: &'a PrimitiveArray<UInt32Type>,
    end: &'a PrimitiveArray<UInt32Type>,
    row_id: &'a PrimitiveArray<Int64Type>,
    i: usize,
}

impl<'a> MergeIterator<'a> {
    pub fn new(recs: &'a RecordBatch) -> MergeIterator<'a> {
        let kind = recs
            .column_by_name("kind")
            .unwrap()
            .as_any()
            .downcast_ref::<DictionaryArray<UInt8Type>>()
            .unwrap();
        let kind_values = kind
            .values()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let chrom_id = recs
            .column_by_name("chrom_id")
            .unwrap()
            .as_any()
            .downcast_ref::<UInt16Array>()
            .unwrap();
        let start = recs
            .column_by_name("start")
            .unwrap()
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        let end = recs
            .column_by_name("end")
            .unwrap()
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        let row_id = recs
            .column_by_name("row_id")
            .unwrap()
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap();
        MergeIterator {
            kind,
            kind_values,
            chrom_id,
            start,
            end,
            row_id,
            i: 0,
        }
    }
}

#[derive(Debug)]
struct Row<'a> {
    kind: &'a str,
    chrom_id: u16,
    start: u32,
    end: u32,
    row_id: i64,
}

impl<'a> Row<'a> {
    pub fn new(kind: &'a str, chrom_id: u16, start: u32, end: u32, row_id: i64) -> Row<'a> {
        Row {
            kind,
            chrom_id,
            start,
            end,
            row_id,
        }
    }
}

impl<'a> HeapItem for Row<'a> {
    type KeyType = u32;

    fn key(&self) -> Self::KeyType {
        self.start
    }
}

impl<'a> Iterator for MergeIterator<'a> {
    type Item = Row<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.i < self.kind.len() {
            let i = self.i;
            self.i += 1;

            let kind = self.kind_values.value(self.kind.key(i).unwrap());
            let chrom_id = self.chrom_id.value(i);
            let start = self.start.value(i);
            let end = self.end.value(i);
            let row_id = self.row_id.value(i);
            Some(Row::new(kind, chrom_id, start, end, row_id))
        } else {
            None
        }
    }
}
