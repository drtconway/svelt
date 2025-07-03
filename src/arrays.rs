use datafusion::arrow::{array::{
    Array, ArrowPrimitiveType, GenericStringBuilder, PrimitiveArray, PrimitiveBuilder, StringArray,
}, datatypes::DataType};

pub trait ArrayWithType<'a>: Array + 'static {
    type ValueType;

    fn value(&'a self, i: usize) -> Self::ValueType;
}

impl<'a, T: ArrowPrimitiveType> ArrayWithType<'a> for PrimitiveArray<T> {
    type ValueType = <T as ArrowPrimitiveType>::Native;

    fn value(&'a self, i: usize) -> Self::ValueType {
        PrimitiveArray::value(self, i)
    }
}

impl<'a> ArrayWithType<'a> for StringArray {
    type ValueType = &'a str;

    fn value(&'a self, i: usize) -> Self::ValueType {
        StringArray::value(self, i)
    }
}

pub trait Builder {
    type ValueType;
    type ArrayType: Array + 'static;

    const DATA_TYPE: DataType;

    fn new() -> Self;

    fn append_value(&mut self, value: Self::ValueType);

    fn append_option(&mut self, value: Option<Self::ValueType>);

    fn finish(&mut self) -> Self::ArrayType;
}

impl<T: ArrowPrimitiveType> Builder for PrimitiveBuilder<T> {
    type ValueType = <T as ArrowPrimitiveType>::Native;
    type ArrayType = PrimitiveArray<T>;

    const DATA_TYPE: DataType = T::DATA_TYPE;

    fn new() -> Self {
        PrimitiveBuilder::new()
    }

    fn append_value(&mut self, value: Self::ValueType) {
        PrimitiveBuilder::append_value(self, value);
    }

    fn append_option(&mut self, value: Option<Self::ValueType>) {
        PrimitiveBuilder::append_option(self, value);
    }

    fn finish(&mut self) -> Self::ArrayType {
        PrimitiveBuilder::finish(self)
    }
}

impl Builder for GenericStringBuilder<i32> {
    type ValueType = String;
    type ArrayType = StringArray;

    const DATA_TYPE: DataType = DataType::Utf8;

    fn new() -> Self {
        GenericStringBuilder::new()
    }

    fn append_value(&mut self, value: Self::ValueType) {
        GenericStringBuilder::append_value(self, value);
    }

    fn append_option(&mut self, value: Option<Self::ValueType>) {
        GenericStringBuilder::append_option(self, value);
    }

    fn finish(&mut self) -> Self::ArrayType {
        GenericStringBuilder::finish(self)
    }
}
