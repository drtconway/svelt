use std::sync::Arc;

use datafusion::arrow::{
    array::RecordBatch,
    datatypes::{Field, Schema, SchemaRef},
};

use crate::arrays::Builder;

pub trait RecordBatchBuilder {
    type RowType;
    type RowOptionType;

    fn append_value(&mut self, row: Self::RowType);

    fn append_option(&mut self, row: Self::RowOptionType);

    fn finish(&mut self) -> RecordBatch;
}

pub struct ColumnBuilder2<Col1Builder: Builder, Col2Builder: Builder> {
    col1: Col1Builder,
    col2: Col2Builder,
    schema: SchemaRef,
}

impl<Col1Builder: Builder, Col2Builder: Builder> ColumnBuilder2<Col1Builder, Col2Builder> {
    pub fn new(col1_name: &str, col1_nullable: bool, col2_name: &str, col2_nullable: bool) -> Self {
        let col1 = Col1Builder::new();
        let col2 = Col2Builder::new();
        let schema = SchemaRef::new(Schema::new(vec![
            Field::new(col1_name, Col1Builder::DATA_TYPE, col1_nullable),
            Field::new(col2_name, Col2Builder::DATA_TYPE, col2_nullable),
        ]));

        ColumnBuilder2 { col1, col2, schema }
    }
}

impl<Col1Builder: Builder, Col2Builder: Builder> RecordBatchBuilder
    for ColumnBuilder2<Col1Builder, Col2Builder>
{
    type RowType = (Col1Builder::ValueType, Col2Builder::ValueType);
    type RowOptionType = (Option<Col1Builder::ValueType>, Option<Col2Builder::ValueType>);

    fn append_value(&mut self, row: Self::RowType) {
        self.col1.append_value(row.0);
        self.col2.append_value(row.1);
    }

    fn append_option(&mut self, row: Self::RowOptionType) {
        self.col1.append_option(row.0);
        self.col2.append_option(row.1);
        
    }
    
    fn finish(&mut self) -> RecordBatch {
        let col1 = self.col1.finish();
        let col2 = self.col2.finish();
        RecordBatch::try_new(self.schema.clone(), vec![Arc::new(col1), Arc::new(col2)]).unwrap()
    }
}
