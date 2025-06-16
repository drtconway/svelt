use std::hash::Hash;
use std::collections::HashMap;

/// An implementation of the disjoint set data structure.
///
/// For more information see https://en.wikipedia.org/wiki/Disjoint-set_data_structure
pub struct DisjointSet<T: Eq + Hash + Copy> {
    parent: HashMap<T, T>,
    rank: HashMap<T, u64>,
}

impl<T: Eq + Hash + Copy> DisjointSet<T> {
    /// Create an empty disjoint set data structure.
    pub fn new() -> DisjointSet<T> {
        return DisjointSet {
            parent: HashMap::new(),
            rank: HashMap::new(),
        };
    }

    /// Find the identifying element of the partition containing `x`.
    pub fn find(&mut self, x: T) -> T {
        if let Some(xp) = self.parent.get(&x) {
            if x != *xp {
                let y = self.find(*xp);
                *self.parent.get_mut(&x).unwrap() = y;
                y
            } else {
                x
            }
        } else {
            self.parent.insert(x, x);
            self.rank.insert(x, 0);
            x
        }
    }

    /// Merge the partition containing `x` with the partition containing `y`.
    pub fn union(&mut self, x: T, y: T) -> T {
        let xr = self.find(x);
        let yr = self.find(y);

        if xr == yr {
            return xr;
        }

        if self.rank.get(&xr).unwrap() < self.rank.get(&yr).unwrap() {
            *self.parent.get_mut(&xr).unwrap() = yr;
            yr
        } else if self.rank.get(&xr).unwrap() > self.rank.get(&yr).unwrap() {
            *self.parent.get_mut(&yr).unwrap() = xr;
            xr
        } else {
            *self.parent.get_mut(&yr).unwrap() = xr;
            *self.rank.get_mut(&xr).unwrap() += 1;
            xr
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn separate_sets() {
        let mut ds = DisjointSet::new();
        assert_eq!(ds.find(4), 4);
        assert_eq!(ds.find(5), 5);
        assert_eq!(ds.find(6), 6);
    }

    #[test]
    fn union_1() {
        let mut ds = DisjointSet::new();
        assert_ne!(ds.find(4), ds.find(5));
        assert_ne!(ds.find(4), ds.find(6));
        ds.union(4, 5);
        assert_eq!(ds.find(4), ds.find(5));
        assert_ne!(ds.find(4), ds.find(6));
    }

    #[test]
    fn union_2() {
        let mut ds = DisjointSet::new();
        assert_ne!(ds.find(4), ds.find(5));
        assert_ne!(ds.find(4), ds.find(6));
        ds.union(4, 5);
        ds.union(5, 6);
        assert_eq!(ds.find(4), ds.find(5));
        assert_eq!(ds.find(4), ds.find(6));
    }
}
