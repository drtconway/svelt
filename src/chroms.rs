use std::collections::HashMap;

pub struct ChromSet {
    names: Vec<String>,
    index: HashMap<String, usize>,
}

impl ChromSet {
    pub fn new() -> Self {
        ChromSet {
            names: Vec::new(),
            index: HashMap::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn add_or_get(&mut self, name: &str) -> usize {
        if let Some(ix) = self.index.get(name) {
            *ix
        } else {
            let n = self.names.len();
            self.names.push(String::from(name));
            self.index.insert(String::from(name), n);
            n
        }
    }

    pub fn name(&self, ix: usize) -> &str {
        &self.names[ix]
    }

    pub fn index(&self, name: &str) -> Option<usize> {
        self.index.get(name).map(|x| *x)
    }
}

impl From<&[&str]> for ChromSet {
    fn from(value: &[&str]) -> Self {
        let names: Vec<String> = value.iter().map(|chrom| String::from(*chrom)).collect();
        let index: HashMap<String, usize> = value
            .iter()
            .enumerate()
            .map(|(ix, chrom)| (String::from(*chrom), ix))
            .collect();
        ChromSet { names, index }
    }
}
