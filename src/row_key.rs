
pub struct RowKey {}

impl RowKey {
    pub fn encode(vix: u32, rn: u32) -> u32 {
        vix + 100 * rn
    }

    pub fn decode(key: u32) -> (u32, u32) {
        (key % 100, key / 100)
    }
}