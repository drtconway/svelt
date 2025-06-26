use std::{
    io::{Error, ErrorKind},
    sync::Arc,
};

use datafusion::{
    arrow::{
        array::{
            DictionaryArray, GenericByteArray, Int32Array, Int64Array, PrimitiveArray,
            PrimitiveBuilder, RecordBatch, StringArray, UInt16Array, UInt32Array,
        },
        datatypes::{
            DataType, Field, GenericStringType, Int32Type, Int64Type, Schema, UInt8Type,
            UInt16Type, UInt32Type,
        },
    },
    common::JoinType,
    prelude::{DataFrame, SessionContext, abs, col, lit},
};

use crate::{
    expressions::prefix_cols,
    heap::{Heap, HeapItem},
    options::MergeOptions,
    row_key::RowKey,
};

pub(super) fn approx_bnd_join(
    orig: DataFrame,
    n: usize,
    options: &MergeOptions,
) -> std::io::Result<DataFrame> {
    let candidates = orig.clone().filter(
        lit(true)
            .and(col("kind").eq(lit("BND")))
            .and(col("vix_count").lt(lit(n as u32))),
    )?;

    let lhs = prefix_cols(candidates.clone(), "lhs")?;
    let rhs = prefix_cols(candidates.clone(), "rhs")?;

    let exact = lhs
        .join(
            rhs,
            JoinType::Inner,
            &["lhs_chrom_id", "lhs_chrom2_id"],
            &["rhs_chrom_id", "rhs_chrom2_id"],
            Some(
                lit(true)
                    .and(col("lhs_row_id").lt(col("rhs_row_id")))
                    .and(col("lhs_row_key").not_eq(col("rhs_row_key")))
                    .and(abs(col("lhs_start") - col("rhs_start")).lt(lit(options.position_window)))
                    .and(abs(col("lhs_end") - col("rhs_end")).lt(lit(options.position_window)))
                    .and(abs(col("lhs_end2") - col("rhs_end2")).lt(lit(options.end2_window)))
                    .and((col("lhs_vix_set") & col("rhs_vix_set")).eq(lit(0))),
            ),
        )?
        .sort(vec![
            (col("lhs_vix_count") + col("rhs_vix_count")).sort(false, false),
            col("lhs_row_key").sort(true, false),
            col("rhs_row_key").sort(true, false),
        ])?
        .select_columns(&["lhs_row_key", "lhs_vix_set", "rhs_row_key", "rhs_vix_set"])?
        .distinct()?;

    Ok(exact)
}

pub(super) async fn approx_near_join(
    orig: DataFrame,
    n: usize,
    options: &MergeOptions,
    ctx: &SessionContext,
) -> std::io::Result<DataFrame> {
    let w = options.position_window as i32;
    let r = options.length_ratio;

    // This is the inner loop test:
    //   1. the starts are within w
    //   2. the ends are within w
    //   3. the length-ratio is >= r
    //   4. the lhs row_id is less than the rhs row_id (avoid symmetric comparisons)
    //   5. the vix (from the row_id) is different - no self-merges.
    let row_is_good = |lhs: &Row<'_>, rhs: &Row<'_>| {
        let lhs_key = RowKey::decode(lhs.row_id as u32);
        let rhs_key = RowKey::decode(rhs.row_id as u32);
        (lhs.start - rhs.start).abs() <= w
            && (lhs.end - rhs.end).abs() <= w
            && lhs.row_id < rhs.row_id
            && lhs.row_key != rhs.row_key
            && lhs_key.0 != rhs_key.0
            && std::cmp::min(lhs.length.abs(), rhs.length.abs()) as f64
                / std::cmp::max(lhs.length.abs(), rhs.length.abs()) as f64
                >= r
    };

    let tbl = orig
        .clone()
        .filter(
            lit(true)
                .and(col("kind").not_eq(lit("BND")))
                .and(col("vix_count").lt(lit(n as u32))),
        )?
        .sort_by(vec![col("kind"), col("chrom_id"), col("start"), col("end")])?;
    let batch = tbl.collect().await?;

    let mut lhs_row_key_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut lhs_vix_set_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut rhs_row_key_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut rhs_vix_set_builder = PrimitiveBuilder::<UInt32Type>::new();

    let mut lhs_heap: Heap<Row<'_>> = Heap::new();
    let mut lhs_itr = batch.iter().flat_map(|recs| MergeIterator::new(recs));
    let mut lhs = lhs_itr.next();

    let mut rhs_heap: Heap<Row<'_>> = Heap::new();
    let mut rhs_itr = batch.iter().flat_map(|recs| MergeIterator::new(recs));
    let mut rhs = rhs_itr.next();

    let mut cached_kind = String::new();
    let mut cached_chrom_id = 0;

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

        if cached_kind != lhs_val.kind {
            cached_kind = String::from(lhs_val.kind);
        }
        if cached_chrom_id != lhs_val.chrom_id {
            cached_chrom_id = lhs_val.chrom_id;
        }

        let lhs_start = lhs_val.start;
        let rhs_start = rhs_val.start;

        if lhs_start <= rhs_start {
            while let Some(rhs_item) = rhs_heap.front() {
                if rhs_item.start + w < lhs_start {
                    rhs_heap.pop();
                } else {
                    break;
                }
            }

            for rhs_item in rhs_heap.iter() {
                assert_eq!(lhs_val.kind, rhs_item.kind);
                assert_eq!(lhs_val.chrom_id, rhs_item.chrom_id);
                assert!(rhs_item.start + w >= lhs_start);
                if !row_is_good(lhs_val, rhs_item) {
                    continue;
                }
                lhs_row_key_builder.append_value(lhs_val.row_key);
                lhs_vix_set_builder.append_value(lhs_val.vix_set);
                rhs_row_key_builder.append_value(rhs_item.row_key);
                rhs_vix_set_builder.append_value(rhs_item.vix_set);
                //log::info!("joining-lhs {:?} and {:?}", lhs_val, rhs_item);
            }

            lhs_heap.push(lhs.take().unwrap());
            lhs = lhs_itr.next();
        } else {
            while let Some(lhs_item) = lhs_heap.front() {
                if lhs_item.start + w < rhs_start {
                    lhs_heap.pop();
                } else {
                    break;
                }
            }

            for lhs_item in lhs_heap.iter() {
                assert_eq!(rhs_val.kind, lhs_item.kind);
                assert_eq!(rhs_val.chrom_id, lhs_item.chrom_id);
                assert!(lhs_item.start + w >= rhs_start);
                if !row_is_good(lhs_item, rhs_val) {
                    continue;
                }
                lhs_row_key_builder.append_value(lhs_item.row_key);
                lhs_vix_set_builder.append_value(lhs_item.vix_set);
                rhs_row_key_builder.append_value(rhs_val.row_key);
                rhs_vix_set_builder.append_value(rhs_val.vix_set);
                //log::info!("joining-rhs {:?} and {:?}", rhs_val, lhs_item);
            }

            rhs_heap.push(rhs.take().unwrap());
            rhs = rhs_itr.next();
        }
    }

    // Pick up the last items on either side.
    // In practice only one side will fire, since
    // we already dropped out of the loop above.

    // Try the LHS first.
    while let Some(lhs_val) = &lhs {
        if lhs_val.kind != cached_kind {
            break;
        }
        if lhs_val.chrom_id != cached_chrom_id {
            break;
        }

        let lhs_start = lhs_val.start;

        while let Some(rhs_item) = rhs_heap.front() {
            if rhs_item.start + w < lhs_start {
                rhs_heap.pop();
            } else {
                break;
            }
        }

        for rhs_item in rhs_heap.iter() {
            assert_eq!(lhs_val.kind, rhs_item.kind);
            assert_eq!(lhs_val.chrom_id, rhs_item.chrom_id);
            assert!(rhs_item.start + w >= lhs_start);
            if !row_is_good(lhs_val, rhs_item) {
                continue;
            }
            lhs_row_key_builder.append_value(lhs_val.row_key);
            lhs_vix_set_builder.append_value(lhs_val.vix_set);
            rhs_row_key_builder.append_value(rhs_item.row_key);
            rhs_vix_set_builder.append_value(rhs_item.vix_set);
            //log::info!("joining-lhs {:?} and {:?}", lhs_val, rhs_item);
        }

        lhs_heap.push(lhs.take().unwrap());
        lhs = lhs_itr.next();
    }

    // Now try the RHS
    while let Some(rhs_val) = &rhs {
        if rhs_val.kind != cached_kind {
            break;
        }
        if rhs_val.chrom_id != cached_chrom_id {
            break;
        }

        let rhs_start = rhs_val.start;

        while let Some(lhs_item) = lhs_heap.front() {
            if lhs_item.start + w < rhs_start {
                lhs_heap.pop();
            } else {
                break;
            }
        }

        for lhs_item in lhs_heap.iter() {
            assert_eq!(rhs_val.kind, lhs_item.kind);
            assert_eq!(rhs_val.chrom_id, lhs_item.chrom_id);
            assert!(lhs_item.start + w >= rhs_start);
            if !row_is_good(lhs_item, rhs_val) {
                continue;
            }
            lhs_row_key_builder.append_value(lhs_item.row_key);
            lhs_vix_set_builder.append_value(lhs_item.vix_set);
            rhs_row_key_builder.append_value(rhs_val.row_key);
            rhs_vix_set_builder.append_value(rhs_val.vix_set);
            //log::info!("joining-rhs {:?} and {:?}", rhs_val, lhs_item);
        }

        rhs_heap.push(rhs.take().unwrap());
        rhs = rhs_itr.next();
    }

    let lhs_row_key_array = lhs_row_key_builder.finish();
    let lhs_vix_set_array = lhs_vix_set_builder.finish();
    let rhs_row_key_array = rhs_row_key_builder.finish();
    let rhs_vix_set_array = rhs_vix_set_builder.finish();

    let schema = Arc::new(Schema::new(vec![
        Field::new("lhs_row_key", DataType::UInt32, false),
        Field::new("lhs_vix_set", DataType::UInt32, false),
        Field::new("rhs_row_key", DataType::UInt32, false),
        Field::new("rhs_vix_set", DataType::UInt32, false),
    ]));

    let recs = RecordBatch::try_new(
        schema,
        vec![
            Arc::new(lhs_row_key_array),
            Arc::new(lhs_vix_set_array),
            Arc::new(rhs_row_key_array),
            Arc::new(rhs_vix_set_array),
        ],
    )
    .map_err(|e| Error::new(ErrorKind::Other, e))?;

    ctx.read_batch(recs)
        .map_err(|e| Error::new(ErrorKind::Other, e))
}

struct MergeIterator<'a> {
    kind: &'a DictionaryArray<UInt8Type>,
    kind_values: &'a GenericByteArray<GenericStringType<i32>>,
    chrom_id: &'a PrimitiveArray<UInt16Type>,
    start: &'a PrimitiveArray<Int32Type>,
    end: &'a PrimitiveArray<Int32Type>,
    length: &'a PrimitiveArray<Int32Type>,
    row_id: &'a PrimitiveArray<Int64Type>,
    row_key: &'a PrimitiveArray<UInt32Type>,
    vix_set: &'a PrimitiveArray<UInt32Type>,
    i: usize,
}

impl<'a> MergeIterator<'a> {
    pub fn new(recs: &'a RecordBatch) -> MergeIterator<'a> {
        let kind = Self::get_array::<DictionaryArray<UInt8Type>>(recs, "kind");
        let kind_values = kind
            .values()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let chrom_id = Self::get_array::<UInt16Array>(recs, "chrom_id");
        let start = Self::get_array::<Int32Array>(recs, "start");
        let end = Self::get_array::<Int32Array>(recs, "end");
        let length = Self::get_array::<Int32Array>(recs, "length");
        let row_id = Self::get_array::<Int64Array>(recs, "row_id");
        let row_key = Self::get_array::<UInt32Array>(recs, "row_key");
        let vix_set = Self::get_array::<UInt32Array>(recs, "vix_set");
        MergeIterator {
            kind,
            kind_values,
            chrom_id,
            start,
            end,
            length,
            row_id,
            row_key,
            vix_set,
            i: 0,
        }
    }

    fn get_array<Type: 'static>(recs: &'a RecordBatch, name: &str) -> &'a Type {
        recs.column_by_name(name)
            .unwrap()
            .as_any()
            .downcast_ref::<Type>()
            .unwrap()
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
            let length = self.length.value(i);
            let row_id = self.row_id.value(i);
            let row_key = self.row_key.value(i);
            let vix_set = self.vix_set.value(i);
            Some(Row::new(
                kind, chrom_id, start, end, length, row_id, row_key, vix_set,
            ))
        } else {
            None
        }
    }
}

#[derive(Debug)]
struct Row<'a> {
    kind: &'a str,
    chrom_id: u16,
    start: i32,
    end: i32,
    length: i32,
    row_id: i64,
    row_key: u32,
    vix_set: u32,
}

impl<'a> Row<'a> {
    pub fn new(
        kind: &'a str,
        chrom_id: u16,
        start: i32,
        end: i32,
        length: i32,
        row_id: i64,
        row_key: u32,
        vix_set: u32,
    ) -> Row<'a> {
        Row {
            kind,
            chrom_id,
            start,
            end,
            length,
            row_id,
            row_key,
            vix_set,
        }
    }
}

impl<'a> HeapItem for Row<'a> {
    type KeyType = i32;

    fn key(&self) -> Self::KeyType {
        self.start
    }
}
