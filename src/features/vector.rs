pub struct MergeVector {
    kmers: Vec<u64>,
    shift: usize,
    index: Vec<usize>,
    postings: Vec<(u32, u32)>,
    toc: Vec<usize>,
}

impl MergeVector {
    pub fn new(k: usize, items: Vec<(u64, Vec<(u32, u32)>)>) -> Self {
        let mut items = items;
        items.sort_by(|lhs, rhs| lhs.0.cmp(&rhs.0));

        let shift = if k > 4 { 2 * k - 8 } else { 0 };
        let mut kmers = Vec::new();
        let mut index = vec![0; 257];
        let mut postings = Vec::new();
        let mut toc = Vec::new();

        toc.push(0);
        for (i, (kmer, mut hits)) in items.into_iter().enumerate() {
            kmers.push(kmer);
            index[1 + (kmer >> shift) as usize] = i;
            postings.append(&mut hits);
            toc.push(postings.len());
        }

        MergeVector {
            kmers,
            shift,
            index,
            postings,
            toc,
        }
    }

    pub fn iter(&self) -> MergeVectorCursor<'_> {
        MergeVectorCursor::new(&self.kmers, self.shift, &self.index, &self.postings, &self.toc)
    }
}

pub struct MergeVectorCursor<'a> {
    kmers: &'a Vec<u64>,
    shift: usize,
    index: &'a Vec<usize>,
    postings: &'a Vec<(u32, u32)>,
    toc: &'a Vec<usize>,
    i: usize,
}

impl<'a> MergeVectorCursor<'a> {
    pub fn new(
        kmers: &'a Vec<u64>,
        shift: usize,
        index: &'a Vec<usize>,
        postings: &'a Vec<(u32, u32)>,
        toc: &'a Vec<usize>,
    ) -> Self {
        MergeVectorCursor {
            kmers,
            shift,
            index,
            postings,
            toc,
            i: 0,
        }
    }

    pub fn seek(&mut self, x: u64) {
        let j = self.index[(x >> self.shift) as usize];
        if j > self.i {
            self.i = j;
        }
        while self.i < self.kmers.len() && self.kmers[self.i] < x {
            self.i += 1;
        }
    }

    pub fn here(&self) -> Option<(u64, &[(u32, u32)])> {
        if self.i < self.kmers.len() {
            let begin = self.toc[self.i];
            let end = self.toc[self.i + 1];
            Some((self.kmers[self.i], &self.postings[begin..end]))
        } else {
            None
        }
    }
}

impl<'a> Iterator for MergeVectorCursor<'a> {
    type Item = (u64, &'a [(u32, u32)]);

    fn next(&mut self) -> Option<Self::Item> {
        if self.i < self.kmers.len() {
            let i = self.i;
            self.i += 1;
            let begin = self.toc[i];
            let end = self.toc[i];
            Some((self.kmers[i], &self.postings[begin..end]))
        } else {
            None
        }
    }
}
