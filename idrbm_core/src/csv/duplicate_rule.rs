/// What to do with duplicate entries?
#[derive(Clone, Copy)]
pub enum DuplicateRule {
    /// The first entry stays.
    FirstWins,
    /// Entries are overwritten, effectively allowing
    /// the last entry within a duplicate group to stay.
    LastWins
}