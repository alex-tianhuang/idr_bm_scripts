use bumpalo::Bump;
use clap::Parser;
use idrbm_core::{csv::{DuplicateRule, read_regions}, datatypes::AMINOACIDS};
use rand::{
    RngExt, SeedableRng,
    distr::slice::Choose,
    rngs::{SmallRng, ThreadRng},
};
use std::path::PathBuf;
/// Program for generating random sequences according
/// to a region CSV file.
#[derive(Parser)]
#[command(verbatim_doc_comment)]
struct Args {
    /// Input regions CSV file.
    input_file: PathBuf,
    /// Output CSV of random sequences.
    output_file: PathBuf,
    /// Number of random sequences to generate.
    #[arg(short)]
    n: usize,
}
fn main() -> anyhow::Result<()> {
    let Args {
        input_file,
        output_file,
        n,
    } = Args::try_parse()?;
    let arena = Bump::new();
    let regions = read_regions(&input_file, DuplicateRule::LastWins, &arena)?;

    let mut writer = csv::Writer::from_path(&output_file)?;
    writer.write_record(&["ProteinID", "RegionID", "VariantID", "Sequence"])?;
    let mut variant_sequence_buffer = String::new();
    let variant_ids = (0..n).map(|i| i.to_string()).collect::<Vec<_>>();

    let rng = SmallRng::from_rng(&mut ThreadRng::default());
    let mut sampler = rng.sample_iter(Choose::new(&AMINOACIDS).unwrap());

    for (protein_id, entries) in regions {
        for &(region_id, region) in entries {
            for variant_id in variant_ids.iter() {
                variant_sequence_buffer.clear();
                variant_sequence_buffer.extend(
                    sampler
                        .by_ref()
                        .take(region.size())
                        .map(|aa| *aa as u8 as char),
                );
                writer.write_record([
                    protein_id,
                    region_id,
                    &variant_id,
                    &variant_sequence_buffer,
                ])?;
            }
        }
    }
    Ok(())
}
