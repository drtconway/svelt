use std::collections::HashMap;


pub struct ChromSet {
    names: Vec<String>,
    index: HashMap<String, usize>
}

impl ChromSet {
    pub fn new() -> Self {
        ChromSet { names: Vec::new(), index: HashMap::new() }
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