
/// A k-length nucleotide sequence represented as a 64-bit integer.
#[derive(PartialEq, Eq, PartialOrd, Ord, Hash, Clone, Debug)]
pub struct Kmer(pub u64);

impl Kmer {
    /// Construct a k-mer from an unsigned integer
    pub fn from_u64(value: u64) -> Kmer {
        Kmer(value)
    }

    /// convert a k-length string into a k-mer.
    ///
    /// If k > 32 or the string contains letters other than
    /// "A", "C", "G", "T", or "U", `None` is returned.
    pub fn make(seq: &str) -> Option<Kmer> {
        if seq.len() > 32 {
            return None;
        }
        let mut x = Kmer(0);
        for c in seq.chars() {
            let b = Kmer::base(c)?;
            x.0 = (x.0 << 2) | b.0;
        }
        Some(x)
    }

    /// Extract k-mers from a string.
    ///
    /// Traverse the string `seq` and put all the valid k-mers in a vector.
    ///
    pub fn make_many(k: usize, seq: &str) -> Vec<Kmer> {
        let mut xs: Vec<Kmer> = Vec::new();
        Self::with_many(k, &seq, |x| xs.push(x.clone()));
        xs
    }

    /// Extract k-mers from a byte sequence and pass them to a closure.
    ///
    pub fn with_many<S, F>(k: usize, seq: &S, mut f: F) -> ()
    where
        S: AsRef<[u8]>,
        F: FnMut(&Kmer),
    {
        let mut i = 0;
        let msk: u64 = (1 << (2 * k)) - 1;
        let mut x = Kmer(0);
        for c in seq.as_ref() {
            match Kmer::byte(*c) {
                None => {
                    i = 0;
                    x.0 = 0;
                }
                Some(b) => {
                    x.0 = (x.0 << 2) | b.0;
                    i += 1;
                    if i == k {
                        x.0 &= msk;
                        f(&x);
                        i -= 1;
                    }
                }
            }
        }
    }

    /// Extract k-mers from both strands of a byte sequence.
    ///
    pub fn with_many_both<S, F>(k: usize, seq: &S, mut f: F) -> ()
    where
        S: AsRef<[u8]>,
        F: FnMut(&Kmer, &Kmer),
    {
        let shift = 2 * (k - 1);
        let msk: u64 = (1 << (2 * k)) - 1;
        let mut x = Kmer(0);
        let mut y = Kmer(0);
        let mut i = 0;
        for c in seq.as_ref() {
            match Kmer::byte(*c) {
                None => {
                    i = 0;
                    x.0 = 0;
                    y.0 = 0;
                }
                Some(b) => {
                    x.0 = (x.0 << 2) | b.0;
                    y.0 = (y.0 >> 2) | ((3 - b.0) << shift);
                    i += 1;
                    if i == k {
                        x.0 &= msk;
                        f(&x, &y);
                        i -= 1;
                    }
                }
            }
        }
    }

    /// Extract k-mers with a zero-based position from both strands of a byte sequence.
    ///
    pub fn with_many_both_pos<S, F>(k: usize, seq: &S, mut f: F) -> ()
    where
        S: AsRef<[u8]>,
        F: FnMut(usize, &Kmer, &Kmer),
    {
        let shift = 2 * (k - 1);
        let msk: u64 = (1 << (2 * k)) - 1;
        let mut x = Kmer(0);
        let mut y = Kmer(0);
        let mut i = 0;
        let mut pos = 0;
        for c in seq.as_ref() {
            match Kmer::byte(*c) {
                None => {
                    i = 0;
                    x.0 = 0;
                    y.0 = 0;
                }
                Some(b) => {
                    x.0 = (x.0 << 2) | b.0;
                    y.0 = (y.0 >> 2) | ((3 - b.0) << shift);
                    i += 1;
                    if i == k {
                        x.0 &= msk;
                        f(pos, &x, &y);
                        i -= 1;
                    }
                }
            }
            pos += 1;
        }
    }
    
    /// Convert a character to a 1-mer.
    #[inline]
    pub fn base(c: char) -> Option<Kmer> {
        match c {
            'A' | 'a' => Some(Kmer(0)),
            'C' | 'c' => Some(Kmer(1)),
            'G' | 'g' => Some(Kmer(2)),
            'T' | 't' | 'U' | 'u' => Some(Kmer(3)),
            _ => None,
        }
    }

    /// Convert an ASCII byte to a 1-mer.
    #[inline]
    pub fn byte(c: u8) -> Option<Kmer> {
        match c {
            b'A' | b'a' => Some(Kmer(0)),
            b'C' | b'c' => Some(Kmer(1)),
            b'G' | b'g' => Some(Kmer(2)),
            b'T' | b't' | b'U' | b'u' => Some(Kmer(3)),
            _ => None,
        }
    }

    /// Create a reversed k-mer.
    #[inline]
    pub fn rev(&self, k: usize) -> Kmer {
        const M2: u64 = 0x3333333333333333;
        const M3: u64 = 0x0F0F0F0F0F0F0F0F;
        const M4: u64 = 0x00FF00FF00FF00FF;
        const M5: u64 = 0x0000FFFF0000FFFF;
        const M6: u64 = 0x00000000FFFFFFFF;

        let mut x = Kmer(self.0);
        x.0 = ((x.0 >> 2) & M2) | ((x.0 & M2) << 2);
        x.0 = ((x.0 >> 4) & M3) | ((x.0 & M3) << 4);
        x.0 = ((x.0 >> 8) & M4) | ((x.0 & M4) << 8);
        x.0 = ((x.0 >> 16) & M5) | ((x.0 & M5) << 16);
        x.0 = ((x.0 >> 32) & M6) | ((x.0 & M6) << 32);
        x.0 >>= 64 - 2 * k;
        x
    }

    /// Create a reverse-complement k-mer.
    #[inline]
    pub fn rev_comp(&self, k: usize) -> Kmer {
        Kmer(!self.0).rev(k)
    }

    /// Convert a k-mer to a string.
    pub fn render(&self, k: usize) -> String {
        let mut s = String::new();
        let mut y = self.rev(k);
        for _i in 0..k {
            match y.0 & 3 {
                0 => {
                    s.push('A');
                }
                1 => {
                    s.push('C');
                }
                2 => {
                    s.push('G');
                }
                3 => {
                    s.push('T');
                }
                _ => {
                    unreachable!();
                }
            }
            y.0 >>= 2;
        }
        s
    }

    /// Compute the Hamming distance between two k-mers.
    ///
    pub fn ham(x: &Kmer, y: &Kmer) -> usize {
        const M1: u64 = 0x5555555555555555;
        let z = x.0 ^ y.0;
        let v = (z | (z >> 1)) & M1;
        v.count_ones() as usize
    }

    /// Compute the frequency of k-mers on the forward strand of a sequence.
    ///
    /// Compute a frequency vector with the number of instances in `seq` of
    /// each of the 4**k possible k-mers.
    ///
    /// Note that this is a dense representation, so the vector will have
    /// 4**k elements.
    pub fn frequency_vector<S>(k: usize, seq: &S) -> Vec<usize>
    where
        S: AsRef<[u8]>,
    {
        let n: usize = 1 << (2 * k);
        let mut v: Vec<usize> = Vec::new();
        v.resize(n, 0);
        Kmer::with_many(k, seq, |x| {
            v[x.0 as usize] += 1;
        });
        v
    }

    /// Compute the frequency of k-mers on both strands of a sequence.
    ///
    /// Compute a frequency vector with the number of instances in `seq` of
    /// each of the 4**k possible k-mers.
    ///
    /// Note that this is a dense representation, so the vector will have
    /// 4**k elements.
    pub fn frequency_vector_both<S>(k: usize, seq: &S) -> Vec<usize>
    where
        S: AsRef<[u8]>,
    {
        let n: usize = 1 << (2 * k);
        let mut v: Vec<usize> = Vec::new();
        v.resize(n, 0);
        Kmer::with_many_both(k, seq, |x, y| {
            v[x.0 as usize] += 1;
            v[y.0 as usize] += 1;
        });
        v
    }
}

/// An iterator over the k-mers drawn from a sequence.
pub struct KmerIterator<'a, Src> 
where
Src: Iterator<Item = &'a u8>
{
    k: usize,
    shift: usize,
    mask: u64,
    x: Kmer,
    y: Kmer,
    i: usize,
    src: Src
}

impl<'a, Src> KmerIterator<'a, Src>
where
Src: Iterator<Item = &'a u8>
{
    /// Create a new iterator yields the k-mers from both strands of a sequence.
    pub fn new(k: usize, src: Src) -> KmerIterator<'a, Src> {
        let shift: usize = 2 * (k - 1);
        let mask: u64 = (1 << (2 * k)) - 1;
        let x: Kmer = Kmer(0);
        let y: Kmer = Kmer(0);
        let i: usize = 0;
        KmerIterator { k, shift, mask, x, y, i, src }
    }
}

impl<'a, Src> Iterator for KmerIterator<'a, Src>
where
Src: Iterator<Item = &'a u8>
{
    type Item = (Kmer, Kmer);
    fn next(&mut self) -> Option<Self::Item> {
        while let Some(c) = self.src.next() {
            match Kmer::byte(*c) {
                None => {
                    self.i = 0;
                    self.x.0 = 0;
                    self.y.0 = 0;
                }
                Some(b) => {
                    self.x.0 = (self.x.0 << 2) | b.0;
                    self.y.0 = (self.y.0 >> 2) | ((3 - b.0) << self.shift);
                    self.i += 1;
                    if self.i == self.k {
                        self.x.0 &= self.mask;
                        self.i -= 1;
                        return Some((self.x.clone(), self.y.clone()));
                    }
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_0a() {
        let seq = "CGAT";
        let ox = Kmer::make(seq);
        assert_eq!(ox, Some(Kmer(0b01100011u64)));
    }

    #[test]
    fn test_0b() {
        let seq = "CGXAT";
        let ox = Kmer::make(seq);
        assert_eq!(ox, None);
    }

    #[test]
    fn test_0c() {
        let seq = "CTTTCTGGGGCTAGAGCAGGCAAACGTGGTACA";
        assert_eq!(seq.len(), 33);
        let ox = Kmer::make(seq);
        assert_eq!(ox, None);
    }

    #[test]
    fn test_1() {
        let seq = "CTTTCTGGGGCTAGAGCAGGCAAACGTGGTACAGTCGACTCCATTCTTTCTTCCTCTGAGACCCCTTCCAGGAATTCAANGGCGCTGGTGAGTCATGAGGCCTCGGAGCAGGGAGTGGTGGTGGTTACATAATTCAGATTAACTCTCAGT";
        let k = 11;
        let xs = Kmer::make_many(k, seq);
        assert_eq!(xs.len(), seq.len() - k + 1 - k);

        let ys = vec![
            "CTTTCTGGGGC",
            "TTTCTGGGGCT",
            "TTCTGGGGCTA",
            "TCTGGGGCTAG",
            "CTGGGGCTAGA",
            "TGGGGCTAGAG",
            "GGGGCTAGAGC",
            "GGGCTAGAGCA",
            "GGCTAGAGCAG",
            "GCTAGAGCAGG",
            "CTAGAGCAGGC",
            "TAGAGCAGGCA",
            "AGAGCAGGCAA",
            "GAGCAGGCAAA",
            "AGCAGGCAAAC",
            "GCAGGCAAACG",
            "CAGGCAAACGT",
            "AGGCAAACGTG",
            "GGCAAACGTGG",
            "GCAAACGTGGT",
            "CAAACGTGGTA",
            "AAACGTGGTAC",
            "AACGTGGTACA",
            "ACGTGGTACAG",
            "CGTGGTACAGT",
            "GTGGTACAGTC",
            "TGGTACAGTCG",
            "GGTACAGTCGA",
            "GTACAGTCGAC",
            "TACAGTCGACT",
            "ACAGTCGACTC",
            "CAGTCGACTCC",
            "AGTCGACTCCA",
            "GTCGACTCCAT",
            "TCGACTCCATT",
            "CGACTCCATTC",
            "GACTCCATTCT",
            "ACTCCATTCTT",
            "CTCCATTCTTT",
            "TCCATTCTTTC",
            "CCATTCTTTCT",
            "CATTCTTTCTT",
            "ATTCTTTCTTC",
            "TTCTTTCTTCC",
            "TCTTTCTTCCT",
            "CTTTCTTCCTC",
            "TTTCTTCCTCT",
            "TTCTTCCTCTG",
            "TCTTCCTCTGA",
            "CTTCCTCTGAG",
            "TTCCTCTGAGA",
            "TCCTCTGAGAC",
            "CCTCTGAGACC",
            "CTCTGAGACCC",
            "TCTGAGACCCC",
            "CTGAGACCCCT",
            "TGAGACCCCTT",
            "GAGACCCCTTC",
            "AGACCCCTTCC",
            "GACCCCTTCCA",
            "ACCCCTTCCAG",
            "CCCCTTCCAGG",
            "CCCTTCCAGGA",
            "CCTTCCAGGAA",
            "CTTCCAGGAAT",
            "TTCCAGGAATT",
            "TCCAGGAATTC",
            "CCAGGAATTCA",
            "CAGGAATTCAA",
            "GGCGCTGGTGA",
            "GCGCTGGTGAG",
            "CGCTGGTGAGT",
            "GCTGGTGAGTC",
            "CTGGTGAGTCA",
            "TGGTGAGTCAT",
            "GGTGAGTCATG",
            "GTGAGTCATGA",
            "TGAGTCATGAG",
            "GAGTCATGAGG",
            "AGTCATGAGGC",
            "GTCATGAGGCC",
            "TCATGAGGCCT",
            "CATGAGGCCTC",
            "ATGAGGCCTCG",
            "TGAGGCCTCGG",
            "GAGGCCTCGGA",
            "AGGCCTCGGAG",
            "GGCCTCGGAGC",
            "GCCTCGGAGCA",
            "CCTCGGAGCAG",
            "CTCGGAGCAGG",
            "TCGGAGCAGGG",
            "CGGAGCAGGGA",
            "GGAGCAGGGAG",
            "GAGCAGGGAGT",
            "AGCAGGGAGTG",
            "GCAGGGAGTGG",
            "CAGGGAGTGGT",
            "AGGGAGTGGTG",
            "GGGAGTGGTGG",
            "GGAGTGGTGGT",
            "GAGTGGTGGTG",
            "AGTGGTGGTGG",
            "GTGGTGGTGGT",
            "TGGTGGTGGTT",
            "GGTGGTGGTTA",
            "GTGGTGGTTAC",
            "TGGTGGTTACA",
            "GGTGGTTACAT",
            "GTGGTTACATA",
            "TGGTTACATAA",
            "GGTTACATAAT",
            "GTTACATAATT",
            "TTACATAATTC",
            "TACATAATTCA",
            "ACATAATTCAG",
            "CATAATTCAGA",
            "ATAATTCAGAT",
            "TAATTCAGATT",
            "AATTCAGATTA",
            "ATTCAGATTAA",
            "TTCAGATTAAC",
            "TCAGATTAACT",
            "CAGATTAACTC",
            "AGATTAACTCT",
            "GATTAACTCTC",
            "ATTAACTCTCA",
            "TTAACTCTCAG",
            "TAACTCTCAGT",
        ];
        assert_eq!(xs.len(), ys.len());
        for i in 0..xs.len() {
            assert_eq!(xs[i].render(k), ys[i]);
        }
    }

    #[test]
    fn test_2() {
        let seq = "CTTTCTGGGGCTAGAGCAGGCAAACGTGGTACAGTCGACTCCATTCTTTCTTCCTCTGAGACCCCTTCCAGGAATTCAANGGCGCTGGTGAGTCATGAGGCCTCGGAGCAGGGAGTGGTGGTGGTTACATAATTCAGATTAACTCTCAGT";
        let k = 11;
        let mut xs = Vec::new();
        Kmer::with_many_both(k, &seq, |x, y| {
            xs.push(x.clone());
            xs.push(y.clone())
        });
        assert_eq!(xs.len(), 2 * (seq.len() - k + 1 - k));
    }

    #[test]
    fn test_3() {
        let seq = "CTTTCTGGGGCTAGAGCAGGCAAACGTGGTACAGTCGACTCCATTCTTTCTTCCTCTGAGACCCCTTCCAGGAATTCAAAGGCGCTGGTGAGTCATGAGGCCTCGGAGCAGGGAGTGGTGGTGGTTACATAATTCAGATTAACTCTCAGT";
        let xs = Kmer::frequency_vector(3, &seq);
        assert_eq!(
            xs,
            vec![
                2, 2, 1, 2, 2, 1, 1, 2, 3, 2, 5, 4, 1, 0, 1, 4, 2, 0, 6, 3, 2, 2, 0, 3, 1, 1, 1, 1,
                1, 5, 3, 4, 1, 2, 6, 1, 3, 1, 1, 2, 3, 4, 3, 5, 1, 2, 5, 1, 2, 2, 1, 0, 4, 3, 2, 5,
                3, 0, 6, 0, 2, 7, 0, 2
            ]
        );
    }

    #[test]
    fn test_4() {
        let seq = "CTTTCTGGGGCTAGAGCAGGCAAACGTGGTACAGTCGACTCCATTCTTTCTTCCTCTGAGACCCCTTCCAGGAATTCAAAGGCGCTGGTGAGTCATGAGGCCTCGGAGCAGGGAGTGGTGGTGGTTACATAATTCAGATTAACTCTCAGT";
        let xs = Kmer::frequency_vector_both(3, &seq);
        assert_eq!(
            xs,
            vec![
                4, 3, 5, 6, 2, 6, 2, 6, 8, 4, 8, 6, 1, 1, 4, 6, 2, 5, 9, 4, 8, 5, 1, 8, 3, 2, 1, 2,
                2, 11, 9, 5, 8, 4, 11, 1, 3, 5, 2, 4, 6, 5, 5, 6, 3, 4, 5, 3, 4, 3, 2, 1, 7, 6, 3,
                8, 7, 3, 8, 2, 4, 8, 2, 4
            ]
        );
    }

    #[test]
    fn test_5() {
        let seq = "CTTTCTGGGGC";
        let k = seq.len();
        let x = Kmer::make(seq).unwrap();
        let y = x.rev_comp(k);
        assert_eq!(y.rev_comp(k), x);
        assert_eq!(y.render(k), "GCCCCAGAAAG");
    }

    #[test]
    fn test_6() {
        let x = Kmer::make("CTTTCTGGGGC").unwrap();
        let y = Kmer::make("CTATGTGGCGC").unwrap();
        assert_eq!(Kmer::ham(&x, &x), 0);
        assert_eq!(Kmer::ham(&y, &y), 0);
        assert_eq!(Kmer::ham(&x, &y), 3);
        assert_eq!(Kmer::ham(&y, &x), 3);
    }

    #[test]
    fn test_7() {
        let seq = "CTTTCTGGGGCTANGGCAAACGTGG";
        let k = 11;
        let mut itr = KmerIterator::new(k, seq.as_bytes().iter());
        let x1 = Kmer::make("CTTTCTGGGGC").unwrap();
        assert_eq!(itr.next(), Some((x1.clone(), x1.rev_comp(k))));
        let x2 = Kmer::make("TTTCTGGGGCT").unwrap();
        assert_eq!(itr.next(), Some((x2.clone(), x2.rev_comp(k))));
        let x3 = Kmer::make("TTCTGGGGCTA").unwrap();
        assert_eq!(itr.next(), Some((x3.clone(), x3.rev_comp(k))));
        let x4 = Kmer::make("GGCAAACGTGG").unwrap();
        assert_eq!(itr.next(), Some((x4.clone(), x4.rev_comp(k))));
        assert_eq!(itr.next(), None);
    }
}
