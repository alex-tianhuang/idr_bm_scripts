use anyhow::Error;
use clap::Parser;
use idrbm_core::datatypes::AMINOACIDS;
use rand::{
    RngExt, SeedableRng,
    distr::slice::Choose,
    rngs::{SmallRng, ThreadRng},
};
use std::path::PathBuf;
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
    let mut writer = csv::Writer::from_path(&output_file)?;
    writer.write_record(&["ProteinID", "RegionID", "VariantID", "Sequence"])?;
    let mut reader = csv::Reader::from_reader(&*contents);
    let mut record = csv::StringRecord::new();
    let rng = SmallRng::from_rng(&mut ThreadRng::default());
    let mut sampler = rng.sample_iter(Choose::new(&AMINOACIDS).unwrap());
    let mut row_buffer = [const { String::new() }; 4];
    let mut variant_ids = (0..n).map(|i| i.to_string()).collect::<Vec<_>>();
    while reader.read_record(&mut record)? {
        // empty stuff
        if record.len() == 0 {
            continue;
        }
        if record.len() != 4 {
            let position = record
                .position()
                .expect("position should be available if the record was successfully read");
            return Err(error(
                &format_args!("expected four parts, got {}", record.len()),
                &contents,
                position,
            ));
        }
        let _: [_; 2] = [0, 1].map(|i| {
            let field = record.get(i).unwrap();
            row_buffer[i].clear();
            row_buffer[i].push_str(field);
            field
        });
        let [start, stop] = [2, 3].map(|i| {
            let field = record.get(i).unwrap();
            field.parse::<usize>().map_err(|e| {
                let position = record
                    .position()
                    .expect("position should be available if the record was successfully read");
                error(&e, &contents, position)
            })
        });
        let start = start?;
        let stop = stop?;
        if start >= stop {
            let position = record
                .position()
                .expect("position should be available if the record was successfully read");
            return Err(error(&"start >= stop", &contents, position))
        }
        for variant_id in variant_ids.iter_mut() {
            let [.., variant_id_buffer, variant_sequence] = &mut row_buffer;
            std::mem::swap(variant_id, variant_id_buffer);
            variant_sequence.clear();
            variant_sequence.extend(
                sampler
                    .by_ref()
                    .take(stop - start)
                    .map(|aa| *aa as u8 as char),
            );
            writer.write_record(&row_buffer)?;
            let [.., variant_id_buffer, _] = &mut row_buffer;
            std::mem::swap(variant_id, variant_id_buffer);
        }
    }
    Ok(())
}
fn error(reason: &dyn std::fmt::Display, contents: &[u8], position: &csv::Position) -> Error {
    let contents = &contents[position.byte() as usize..];
    let next_newline = contents
        .iter()
        .position(|b| *b == b'\n')
        .unwrap_or(contents.len());
    let line = &contents[..next_newline];
    Error::msg(format!(
        "unable to read line {} ({}): {}",
        position.line(),
        reason,
        String::from_utf8_lossy(line)
    ))
}
