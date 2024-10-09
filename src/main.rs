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

    #[arg(
		long,
		default_value_t = 1,
		value_parser = clap::value_parser!(u32).range(1..),
	)]
    page: u32,

    #[arg(
	    long,
		default_value="1:500,1:1000",
		num_args=1..,
		value_delimiter=',',
		value_parser=parse_scale,
	)]
    scales: Vec<u32>,
}

pub fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    let matcher = Matcher::new(args.points, args.symbols)?;
    let dpi = 300; // TODO: Retrieve image resulution from TIFF tags.
    if let Some(tr) = matcher.find_transform(dpi, &args.scales) {
        write_geotiff(args.image, args.page, &tr, args.output)?;
    }

    Ok(())
}

// Parses a single value of the `--scales` command line argument.
//
// # Examples
//
// ```
// assert_eq!(parse_scales("1:200").unwrap(), 200);
// assert_eq!(parse_scales("1:blah").is_err(), true);
// ```
fn parse_scale(s: &str) -> Result<u32, String> {
    if let Some(stripped) = s.strip_prefix("1:") {
        if let Ok(num) = stripped.parse::<u32>() {
            return Ok(num);
        }
    }
    Err("scale must be formed like \"1:500\" or similar")?
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn can_parse_scale() {
        assert_eq!(parse_scale("1:200").unwrap(), 200);
        assert_eq!(parse_scale("").is_err(), true);
        assert_eq!(parse_scale("1:").is_err(), true);
        assert_eq!(parse_scale("blah").is_err(), true);
    }
}
