pub trait HeapItem {
    type KeyType: PartialEq + PartialOrd;

    fn key(&self) -> Self::KeyType;
}

pub struct Heap<Item: HeapItem> {
    items: Vec<Item>,
}

impl<Item: HeapItem> Heap<Item> {
    pub fn new() -> Heap<Item> {
        Heap { items: Vec::new() }
    }

    pub fn clear(&mut self) {
        self.items.clear();
    }

    pub fn push(&mut self, item: Item) {
        self.items.push(item);
        let n = self.items.len();
        self.upheap(n);
    }

    pub fn pop(&mut self) -> Option<Item> {
        let n = self.items.len();
        if n > 0 {
            if n > 1 {
                self.items.swap(0, n - 1);
            }
            let res = self.items.pop();
            self.downheap(1);
            res
        } else {
            None
        }
    }

    pub fn front(&self) -> Option<&Item> {
        self.items.first()
    }

    pub fn iter(&self) -> impl Iterator<Item = &Item> {
        self.items.iter()
    }

    pub fn change_front(&mut self, item: Item) {
        self.items[0] = item;
        self.downheap(1);
    }

    fn heapify(&mut self) {
        let n = self.items.len() + 1;
        for i in 1..n {
            self.upheap(i);
        }
    }

    fn upheap(&mut self, i: usize) {
        let mut i = i;
        let mut p = i / 2;
        while p >= 1 {
            let i_key = self.items[i - 1].key();
            let p_key = self.items[p - 1].key();
            if i_key < p_key {
                self.items.swap(i - 1, p - 1);
                i = p;
                p = i / 2;
            } else {
                break;
            }
        }
    }

    fn downheap(&mut self, p: usize) {
        let n = self.items.len();
        let mut p = p;
        let mut c = 2 * p;
        while c <= n {
            let mut cm1_key = self.items[c - 1].key();
            if c < n {
                let c_key = self.items[c].key();
                if c_key < cm1_key {
                    c += 1;
                    cm1_key = c_key;
                }
            }
            let p_key = self.items[p - 1].key();
            if cm1_key < p_key {
                self.items.swap(c - 1, p - 1);
                p = c;
                c = 2 * p
            } else {
                break;
            }
        }
    }
}

impl<Item: HeapItem> From<Vec<Item>> for Heap<Item> {
    fn from(value: Vec<Item>) -> Self {
        let mut heap = Heap { items: value };
        heap.heapify();
        heap
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug, PartialEq, PartialOrd)]
    struct X(i64);

    impl HeapItem for X {
        type KeyType = i64;

        fn key(&self) -> Self::KeyType {
            self.0
        }
    }

    #[test]
    fn test_heap_1() {
        let xs: Vec<i64> = vec![4, 7, 11, 2, 5];
        let items: Vec<X> = xs.iter().map(|x| X(*x)).collect();
        let mut h = Heap::from(items);
        assert_eq!(h.front(), Some(&X(2)));
        h.change_front(X(15));
        assert_eq!(h.front(), Some(&X(4)));
    }

    #[test]
    fn test_heap_2() {
        let xs: Vec<i64> = vec![4, 7, 11, 2, 5];
        let items: Vec<X> = xs.iter().map(|x| X(*x)).collect();
        let mut h = Heap::from(items);
        assert_eq!(h.pop(), Some(X(2)));
        assert_eq!(h.front(), Some(&X(4)));
        h.push(X(3));
        assert_eq!(h.front(), Some(&X(3)));
    }
}
