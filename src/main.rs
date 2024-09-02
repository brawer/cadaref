use clap::Parser;
use std::error::Error;
use std::path::PathBuf;

use cadaref::{write_geotiff, Matcher};

#[derive(Parser, Debug)]
struct Args {
    #[arg(long, value_hint=clap::ValueHint::FilePath)]
    points: PathBuf,

    #[arg(long, value_hint=clap::ValueHint::FilePath)]
    symbols: PathBuf,

    #[arg(short, long, value_hint=clap::ValueHint::FilePath)]
    output: PathBuf,

    #[arg(value_hint=clap::ValueHint::FilePath)]
    image: PathBuf,
}

pub fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    let matcher = Matcher::new(args.points, args.symbols)?;

    let dpi = 300.0;
    let pixels_per_meter = dpi / 2.54 * 100.0;
    let meters_per_pixel = 500.0 / pixels_per_meter;
    if let Some((tr, score)) = matcher.find_matches(meters_per_pixel) {
        println!("{:?} {:?}", score, tr);
        write_geotiff(args.image, &tr, args.output)?;
    }

    Ok(())
}
