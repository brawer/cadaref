use clap::Parser;
use std::error::Error;
use std::path::PathBuf;

use cadaref::Matcher;

#[derive(Parser, Debug)]
struct Args {
    #[arg(long, value_hint=clap::ValueHint::FilePath)]
    points: PathBuf,

    #[arg(long, value_hint=clap::ValueHint::FilePath)]
    symbols: PathBuf,
}

pub fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    let matcher = Matcher::new(args.points, args.symbols)?;

    let dpi = 300.0;
    let pixels_per_meter = dpi / 2.54 * 100.0;
    let meters_per_pixel = 500.0 / pixels_per_meter;
    println!("Hello world {}", meters_per_pixel);

    if let Some((t, score)) = matcher.find_matches(meters_per_pixel) {
        println!("{:?} {:?}", score, t);
    }

    Ok(())
}
