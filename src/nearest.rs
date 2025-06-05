use crate::errors::SveltError;

pub struct Neardex<Value> {
    items: Vec<(u32, Value)>,
}

impl<Value> Neardex<Value> {
    pub fn new(items: Vec<(u32, Value)>) -> std::result::Result<Neardex<Value>, SveltError> {
        let mut res = Neardex { items };
        res.items.sort_by(|lhs, rhs| lhs.0.cmp(&rhs.0));
        for i in 1..res.items.len() {
            if res.items[i - 1].0 == res.items[i].0 {
                return Err(SveltError::NeardexDuplicate(res.items[i].0));
            }
        }
        Ok(res)
    }

    pub fn len(&self) -> usize {
        self.items.len()
    }

    pub fn nearest(&self, x: u32) -> Option<&(u32, Value)> {
        let i = self.lower_bound(x);
        let n = self.len();
        match (i > 0, i < n) {
            (true, true) => {
                let x0 = self.items[i - 1].0;
                let x1 = self.items[i].0;
                if x - x0 <= x1 - x {
                    Some(&self.items[i - 1])
                } else {
                    Some(&self.items[i])
                }
            }
            (true, false) => Some(&self.items[i - 1]),
            (false, true) => Some(&self.items[i]),
            (false, false) => None,
        }
    }

    pub fn within(&self, x: u32, d: u32) -> Vec<&(u32, Value)> {
        let mut res = Vec::new();
        let i = self.lower_bound(x);
        let mut down = self.items[..i]
            .iter()
            .rev()
            .map(|item| (x - item.0, item))
            .take_while(|item| item.0 < d);
        let mut up = self.items[i..]
            .iter()
            .map(|item| (item.0 - x, item))
            .take_while(|item| item.0 < d);

        let mut next_down = down.next();
        let mut next_up = up.next();
        while let (Some(lhs), Some(rhs)) = (&next_down, &next_up) {
            if lhs.0 <= rhs.0 {
                res.push(lhs.1);
                next_down = down.next();
            } else {
                res.push(rhs.1);
                next_up = up.next();
            }
        }
        while let Some(lhs) = &next_down {
            res.push(lhs.1);
            next_down = down.next();
        }
        while let Some(rhs) = &next_up {
            res.push(rhs.1);
            next_up = up.next();
        }
        res
    }

    fn lower_bound(&self, x: u32) -> usize {
        let mut first = 0;
        let mut count = self.items.len();
        while count > 0 {
            let step = count / 2;
            let i = first + step;
            if self.items[i].0 < x {
                first = i + 1;
                count -= step + 1;
            } else {
                count = step;
            }
        }
        first
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nearest() {
        let items = vec![(2, 'a'), (3, 'b'), (5, 'c'), (8, 'd'), (11, 'e')];
        let dex = Neardex::new(items).unwrap();
        assert_eq!(dex.nearest(2), Some(&(2, 'a')));
        assert_eq!(dex.nearest(4), Some(&(3, 'b')));
        assert_eq!(dex.nearest(7), Some(&(8, 'd')));
    }

    #[test]
    fn test_within() {
        let items = vec![(2, 'a'), (3, 'b'), (5, 'c'), (8, 'd'), (11, 'e')];
        let dex = Neardex::new(items).unwrap();
        assert_eq!(dex.within(4, 2), vec![&(3, 'b'), &(5, 'c')]);
        assert_eq!(dex.within(4, 3), vec![&(3, 'b'), &(5, 'c'), &(2, 'a')]);
    }
}
