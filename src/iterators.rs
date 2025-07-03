use std::iter::{Zip, zip};

use datafusion::arrow::array::{Array, RecordBatch};

use crate::arrays::ArrayWithType;

pub struct ColumnIterator2<'a, Col1Type: ArrayWithType<'a>, Col2Type: ArrayWithType<'a>> {
    col1: &'a Col1Type,
    col2: &'a Col2Type,
    i: usize,
}

impl<'a, Col1Type: ArrayWithType<'a>, Col2Type: ArrayWithType<'a>>
    ColumnIterator2<'a, Col1Type, Col2Type>
{
    pub fn new(recs: &'a RecordBatch, col1_name: &str, col2_name: &str) -> Self {
        let col1 = get_array::<Col1Type>(recs, col1_name);
        let col2 = get_array::<Col2Type>(recs, col2_name);

        ColumnIterator2 { col1, col2, i: 0 }
    }
}

impl<'a, Col1Type: ArrayWithType<'a>, Col2Type: ArrayWithType<'a>> Iterator
    for ColumnIterator2<'a, Col1Type, Col2Type>
{
    type Item = (Col1Type::ValueType, Col2Type::ValueType);

    fn next(&mut self) -> Option<Self::Item> {
        if self.i < self.col1.len() {
            let value = (self.col1.value(self.i), self.col2.value(self.i));
            self.i += 1;
            Some(value)
        } else {
            None
        }
    }
}

pub struct ColumnIterator3<
    'a,
    Col1Type: ArrayWithType<'a>,
    Col2Type: ArrayWithType<'a>,
    Col3Type: ArrayWithType<'a>,
> {
    col1: &'a Col1Type,
    col2: &'a Col2Type,
    col3: &'a Col3Type,
    i: usize,
}

impl<'a, Col1Type: ArrayWithType<'a>, Col2Type: ArrayWithType<'a>, Col3Type: ArrayWithType<'a>>
    ColumnIterator3<'a, Col1Type, Col2Type, Col3Type>
{
    pub fn new(recs: &'a RecordBatch, col1_name: &str, col2_name: &str, col3_name: &str) -> Self {
        let col1 = get_array::<Col1Type>(recs, col1_name);
        let col2 = get_array::<Col2Type>(recs, col2_name);
        let col3 = get_array::<Col3Type>(recs, col3_name);

        ColumnIterator3 {
            col1,
            col2,
            col3,
            i: 0,
        }
    }
}

impl<'a, Col1Type: ArrayWithType<'a>, Col2Type: ArrayWithType<'a>, Col3Type: ArrayWithType<'a>>
    Iterator for ColumnIterator3<'a, Col1Type, Col2Type, Col3Type>
{
    type Item = (
        Col1Type::ValueType,
        Col2Type::ValueType,
        Col3Type::ValueType,
    );

    fn next(&mut self) -> Option<Self::Item> {
        if self.i < self.col1.len() {
            let value = (
                self.col1.value(self.i),
                self.col2.value(self.i),
                self.col3.value(self.i),
            );
            self.i += 1;
            Some(value)
        } else {
            None
        }
    }
}

pub struct ColumnIterator4<
    'a,
    Col1Type: ArrayWithType<'a>,
    Col2Type: ArrayWithType<'a>,
    Col3Type: ArrayWithType<'a>,
    Col4Type: ArrayWithType<'a>,
> {
    iter: Zip<ColumnIterator2<'a, Col1Type, Col2Type>, ColumnIterator2<'a, Col3Type, Col4Type>>,
}

impl<
    'a,
    Col1Type: ArrayWithType<'a>,
    Col2Type: ArrayWithType<'a>,
    Col3Type: ArrayWithType<'a>,
    Col4Type: ArrayWithType<'a>,
> ColumnIterator4<'a, Col1Type, Col2Type, Col3Type, Col4Type>
{
    pub fn new(
        recs: &'a RecordBatch,
        col1_name: &str,
        col2_name: &str,
        col3_name: &str,
        col4_name: &str,
    ) -> Self {
        ColumnIterator4 {
            iter: zip(
                ColumnIterator2::new(recs, col1_name, col2_name),
                ColumnIterator2::new(recs, col3_name, col4_name),
            ),
        }
    }
}

impl<
    'a,
    Col1Type: ArrayWithType<'a>,
    Col2Type: ArrayWithType<'a>,
    Col3Type: ArrayWithType<'a>,
    Col4Type: ArrayWithType<'a>,
> Iterator for ColumnIterator4<'a, Col1Type, Col2Type, Col3Type, Col4Type>
{
    type Item = (
        (Col1Type::ValueType, Col2Type::ValueType),
        (Col3Type::ValueType, Col4Type::ValueType),
    );

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

pub struct ColumnIterator5<
    'a,
    Col1Type: ArrayWithType<'a>,
    Col2Type: ArrayWithType<'a>,
    Col3Type: ArrayWithType<'a>,
    Col4Type: ArrayWithType<'a>,
    Col5Type: ArrayWithType<'a>,
> {
    iter: Zip<
        ColumnIterator3<'a, Col1Type, Col2Type, Col3Type>,
        ColumnIterator2<'a, Col4Type, Col5Type>,
    >,
}

impl<
    'a,
    Col1Type: ArrayWithType<'a>,
    Col2Type: ArrayWithType<'a>,
    Col3Type: ArrayWithType<'a>,
    Col4Type: ArrayWithType<'a>,
    Col5Type: ArrayWithType<'a>,
> ColumnIterator5<'a, Col1Type, Col2Type, Col3Type, Col4Type, Col5Type>
{
    pub fn new(
        recs: &'a RecordBatch,
        col1_name: &str,
        col2_name: &str,
        col3_name: &str,
        col4_name: &str,
        col5_name: &str,
    ) -> Self {
        ColumnIterator5 {
            iter: zip(
                ColumnIterator3::new(recs, col1_name, col2_name, col3_name),
                ColumnIterator2::new(recs, col4_name, col5_name),
            ),
        }
    }
}

impl<
    'a,
    Col1Type: ArrayWithType<'a>,
    Col2Type: ArrayWithType<'a>,
    Col3Type: ArrayWithType<'a>,
    Col4Type: ArrayWithType<'a>,
    Col5Type: ArrayWithType<'a>,
> Iterator for ColumnIterator5<'a, Col1Type, Col2Type, Col3Type, Col4Type, Col5Type>
{
    type Item = (
        (
            Col1Type::ValueType,
            Col2Type::ValueType,
            Col3Type::ValueType,
        ),
        (Col4Type::ValueType, Col5Type::ValueType),
    );

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

fn get_array<'a, Type: 'static>(recs: &'a RecordBatch, name: &str) -> &'a Type {
    if false {
        log::info!("getting {}", name);
    }
    recs.column_by_name(name)
        .unwrap()
        .as_any()
        .downcast_ref::<Type>()
        .unwrap()
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use datafusion::arrow::array::{
        Int16Array, RecordBatch, StringArray, UInt16Array, UInt32Array, UInt64Array,
    };
    use datafusion::arrow::datatypes::{DataType, Field, Schema};

    use super::*;

    #[tokio::test]
    async fn test_col2_1() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("a", DataType::UInt32, false),
            Field::new("b", DataType::Int16, false),
        ]));

        let recs = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt32Array::from(vec![1, 2, 3, 4])),
                Arc::new(Int16Array::from(vec![-1, -2, -3, -4])),
            ],
        )
        .unwrap();

        let mut itr: ColumnIterator2<'_, UInt32Array, Int16Array> =
            ColumnIterator2::new(&recs, "a", "b");
        assert_eq!(itr.next(), Some((1, -1)));
        assert_eq!(itr.next(), Some((2, -2)));
        assert_eq!(itr.next(), Some((3, -3)));
        assert_eq!(itr.next(), Some((4, -4)));
        assert_eq!(itr.next(), None);
    }

    #[tokio::test]
    async fn test_col3_1() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("a", DataType::UInt32, false),
            Field::new("b", DataType::Int16, false),
            Field::new("c", DataType::Utf8, false),
        ]));

        let recs = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt32Array::from(vec![1, 2, 3, 4])),
                Arc::new(Int16Array::from(vec![-1, -2, -3, -4])),
                Arc::new(StringArray::from(vec![
                    String::from("i"),
                    String::from("ii"),
                    String::from("iii"),
                    String::from("iv"),
                ])),
            ],
        )
        .unwrap();

        let mut itr: ColumnIterator3<'_, UInt32Array, Int16Array, StringArray> =
            ColumnIterator3::new(&recs, "a", "b", "c");
        assert_eq!(itr.next(), Some((1, -1, "i")));
        assert_eq!(itr.next(), Some((2, -2, "ii")));
        assert_eq!(itr.next(), Some((3, -3, "iii")));
        assert_eq!(itr.next(), Some((4, -4, "iv")));
        assert_eq!(itr.next(), None);
    }

    #[tokio::test]
    async fn test_col5_1() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("a", DataType::UInt32, false),
            Field::new("b", DataType::Int16, false),
            Field::new("c", DataType::Utf8, false),
            Field::new("d", DataType::UInt16, false),
            Field::new("e", DataType::UInt64, false),
        ]));

        let recs = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt32Array::from(vec![1, 2, 3, 4])),
                Arc::new(Int16Array::from(vec![-1, -2, -3, -4])),
                Arc::new(StringArray::from(vec![
                    String::from("i"),
                    String::from("ii"),
                    String::from("iii"),
                    String::from("iv"),
                ])),
                Arc::new(UInt16Array::from(vec![1 << 1, 1 << 2, 1 << 3, 1 << 4])),
                Arc::new(UInt64Array::from(vec![1 << 8, 1 << 9, 1 << 10, 1 << 11])),
            ],
        )
        .unwrap();

        let mut itr: ColumnIterator5<
            '_,
            UInt32Array,
            Int16Array,
            StringArray,
            UInt16Array,
            UInt64Array,
        > = ColumnIterator5::new(&recs, "a", "b", "c", "d", "e");

        assert_eq!(itr.next(), Some(((1, -1, "i"), (2, 256))));
        assert_eq!(itr.next(), Some(((2, -2, "ii"), (4, 512))));
        assert_eq!(itr.next(), Some(((3, -3, "iii"), (8, 1024))));
        assert_eq!(itr.next(), Some(((4, -4, "iv"), (16, 2048))));
        assert_eq!(itr.next(), None);
    }
}
