use clap::Parser;
use idrbm_core::{datatypes::AMINOACIDS, utils::ParseRegions};
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
    regions: PathBuf,
    /// Output CSV of random sequences.
    output_file: PathBuf,
    /// Number of random sequences to generate.
    #[arg(short)]
    n: usize,
}
fn main() -> anyhow::Result<()> {
    let Args {
        regions,
        output_file,
        n,
    } = Args::try_parse()?;
    let contents = std::fs::read(&regions)?;
    let mut reader = ParseRegions::new(&contents);
    let mut record = csv::StringRecord::new();

    let mut writer = csv::Writer::from_path(&output_file)?;
    writer.write_record(&["ProteinID", "RegionID", "VariantID", "Sequence"])?;
    let mut row_buffer = [const { String::new() }; 4];
    let mut variant_ids = (0..n).map(|i| i.to_string()).collect::<Vec<_>>();

    let rng = SmallRng::from_rng(&mut ThreadRng::default());
    let mut sampler = rng.sample_iter(Choose::new(&AMINOACIDS).unwrap());

    while let Some(notification) = reader.next(&mut record) {
        let record = notification?;
        for (buf, src) in row_buffer
            .iter_mut()
            .zip([record.protein_id, record.region_id])
        {
            buf.clear();
            buf.push_str(src);
        }
        for variant_id in variant_ids.iter_mut() {
            let [.., variant_id_buffer, variant_sequence] = &mut row_buffer;
            std::mem::swap(variant_id, variant_id_buffer);
            variant_sequence.clear();
            variant_sequence.extend(
                sampler
                    .by_ref()
                    .take(record.region_len())
                    .map(|aa| *aa as u8 as char),
            );
            writer.write_record(&row_buffer)?;
            let [.., variant_id_buffer, _] = &mut row_buffer;
            std::mem::swap(variant_id, variant_id_buffer);
        }
    }
    Ok(())
}
