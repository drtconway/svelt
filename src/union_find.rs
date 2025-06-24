use std::sync::{
    Arc,
    atomic::{AtomicU64, Ordering},
};

#[derive(Clone)]
pub struct UnionFind {
    data: Arc<Vec<AtomicU64>>,
}

impl UnionFind {
    pub fn new(n: usize) -> UnionFind {
        let data = Arc::new(
            (0..n)
                .into_iter()
                .map(|i| AtomicU64::new(i as u64))
                .collect::<Vec<AtomicU64>>(),
        );
        UnionFind { data }
    }

    pub fn find(&self, id: u32) -> u32 {
        let mut id = id;
        while id != self.parent(id) {
            let x = self.data[id as usize].load(Ordering::Relaxed);
            let x_id = (x & 0xFFFFFFFF) as u32;
            let y_id = self.parent(x_id);
            let x_rank = x & 0xFFFFFFFF00000000;
            let y = x_rank | (y_id as u64);
            if y != x {
                let _ = self.data[id as usize].compare_exchange_weak(
                    x,
                    y,
                    Ordering::SeqCst,
                    Ordering::Relaxed,
                );
            }
            id = y_id;
        }
        id
    }

    pub fn union(&self, lhs_id: u32, rhs_id: u32) -> u32 {
        let mut lhs_id = lhs_id;
        let mut rhs_id = rhs_id;
        loop {
            lhs_id = self.find(lhs_id);
            rhs_id = self.find(rhs_id);

            if lhs_id == rhs_id {
                return lhs_id;
            }

            let mut lhs_rank = self.rank(lhs_id);
            let mut rhs_rank = self.rank(rhs_id);

            if lhs_rank > rhs_rank || (lhs_rank == rhs_rank && lhs_id < rhs_id) {
                std::mem::swap(&mut lhs_id, &mut rhs_id);
                std::mem::swap(&mut lhs_rank, &mut rhs_rank);
            }

            let old = ((lhs_rank as u64) << 32) | (lhs_id as u64);
            let new = ((lhs_rank as u64) << 32) | (rhs_id as u64);

            if let Err(_) = self.data[lhs_id as usize].compare_exchange(
                old,
                new,
                Ordering::SeqCst,
                Ordering::Relaxed,
            ) {
                continue;
            }

            if lhs_rank == rhs_rank {
                let old = ((rhs_rank as u64) << 32) | (rhs_id as u64);
                let new = (((rhs_rank + 1) as u64) << 32) | (rhs_id as u64);
                let _ = self.data[rhs_id as usize].compare_exchange_weak(
                    old,
                    new,
                    Ordering::SeqCst,
                    Ordering::Relaxed,
                );
            }

            break;
        }

        rhs_id
    }

    fn parent(&self, id: u32) -> u32 {
        (self.data[id as usize].load(Ordering::Relaxed) & 0xFFFFFFFF) as u32
    }

    fn rank(&self, id: u32) -> u32 {
        (self.data[id as usize].load(Ordering::Relaxed) >> 32 & 0xFFFFFFFF) as u32
    }
}

#[cfg(test)]
mod tests {
    use rand::{SeedableRng, rngs::StdRng, seq::SliceRandom};

    use crate::disjoint_set::DisjointSet;

    use super::*;

    #[test]
    fn basic_union_find() {
        let s = 19;
        let mut rng = StdRng::seed_from_u64(s);
        let n = 100;
        let mut xs: Vec<u32> = (0..n).collect();
        xs.shuffle(&mut rng);
        let pairs1: Vec<(u32, u32)> = xs.chunks(2).map(|c| (c[0], c[1])).collect();

        let mut ds = DisjointSet::new();
        let uf = UnionFind::new(n as usize);

        for (x, y) in pairs1.iter() {
            let x = *x;
            let y = *y;

            assert_eq!(uf.find(x), x);
            assert_eq!(uf.find(y), y);
        }

        for (x, y) in pairs1.iter() {
            let x = *x;
            let y = *y;

            ds.union(x, y);
            uf.union(x, y);
        }

        for (x, y) in pairs1.iter() {
            let x = *x;
            let y = *y;

            assert_eq!(ds.find(x), ds.find(y));
            assert_eq!(uf.find(x), uf.find(y));
        }
    }
}
