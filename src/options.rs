pub struct Options {
    /// Allowed distance for merging of events
    pub position_window: u32,

    /// Allowed end2 distance for merging BND events
    pub end2_window: u32,

    /// Minimum length ratio for merging two events
    pub length_ratio: f64,
}

impl Options {
    pub fn set_position_window(mut self, position_window: u32) -> Options {
        self.position_window = position_window;
        self
    }

    pub fn set_end2_window(mut self, end2_window: u32) -> Options {
        self.end2_window = end2_window;
        self
    }

    pub fn set_length_ratio(mut self, length_ratio: f64) -> Options {
        self.length_ratio = length_ratio;
        self
    }
}

impl Default for Options {
    fn default() -> Self {
        Self {
            position_window: 25,
            end2_window: 150,
            length_ratio: 0.9,
        }
    }
}
