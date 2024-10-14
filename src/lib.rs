use bisection::bisect_left_by;
use csv::Reader;
use gdal::{raster::RasterCreationOptions, Dataset, DriverManager, GeoTransform, Metadata};
use gdal_sys::{GDALGCPsToGeoTransform, GDAL_GCP};
use rstar::{primitives::GeomWithData, RTree};
use serde::Deserialize;
use std::error::Error;
use std::{fs::File, path::PathBuf};

#[derive(Debug, Deserialize, PartialEq, PartialOrd)]
struct Symbol {
    x: f64,
    y: f64,
}

#[derive(Debug, Deserialize, PartialEq, PartialOrd)]
struct Point {
    id: String,
    x: f64,
    y: f64,
}

type PointNode = GeomWithData<[f64; 2], usize>;

type PointDistances = Vec<(i32, u32, u32)>; // millimeters, p1, p2

pub struct Matcher {
    points: Vec<Point>,
    point_distances: PointDistances,
    point_tree: RTree<PointNode>,
    symbols: Vec<Symbol>,
}

#[derive(Copy, Clone, Debug)]
pub struct Score {
    num_matches: usize,
    avg_dist: f64,
}

impl Score {
    fn better_than(&self, other: &Self) -> bool {
        if self.num_matches > other.num_matches {
            return true;
        }
        if self.num_matches < other.num_matches {
            return false;
        }
        self.avg_dist < other.avg_dist
    }
}

impl Matcher {
    pub fn new(points: PathBuf, symbols: PathBuf) -> Result<Self, Box<dyn Error>> {
        let points = read_points(points)?;
        let symbols = read_symbols(symbols)?;

        let nodes = points
            .iter()
            .enumerate()
            .map(|(i, p)| PointNode::new([p.x, p.y], i))
            .collect();
        let point_tree = RTree::bulk_load(nodes);
        let point_distances = build_point_distances(&points);

        let matcher = Matcher {
            points,
            point_distances,
            point_tree,
            symbols,
        };
        Ok(matcher)
    }

    // Tolerance for considering a potential match as inlier. If a point
    // is <= 1 meter away from the projected location, it is considered as
    // an inlier.
    const MAX_DISTANCE_MM: i32 = 1000;

    pub fn find_transform(&self, dpi: f64, scales: &Vec<u32>) -> Option<GeoTransform> {
        let mut best_transform = [0.0; 6];
        let mut best_score = Score {
            num_matches: 0,
            avg_dist: 0.0,
        };

        for scale in scales {
            if let Some((tr, score)) = self.find(dpi, *scale) {
                if score.better_than(&best_score) {
                    best_transform = tr;
                    best_score = score;
                }
            }
        }

        if best_score.num_matches > 0 {
            Some(best_transform)
        } else {
            None
        }
    }

    fn find(&self, dpi: f64, scale: u32) -> Option<(GeoTransform, Score)> {
        let pixels_per_meter = dpi / 2.54 * 100.0;
        let meters_per_pixel = (scale as f64) / pixels_per_meter;

        let num_symbols = self.symbols.len();

        let mut best_transform = [0.0; 6];
        let mut best_score = Score {
            num_matches: 0,
            avg_dist: 0.0,
        };

        const MAX_MATCHES: usize = 30;
        for i in 0..num_symbols.min(MAX_MATCHES) {
            let mut start = i + 1;
            if num_symbols - start > MAX_MATCHES {
                start = num_symbols - MAX_MATCHES;
            }

            for j in start..num_symbols {
                let s1 = &self.symbols[i];
                let s2 = &self.symbols[j];
                let dist_mm = s1.distance_mm(s2, meters_per_pixel);

                // If the two map symbols S1 and S2 are very near each other,
                // we ignore the pair because the computated transform would
                // not be accurate.
                if dist_mm < 10000 {
                    // 10 meters
                    continue;
                }

                let min_dist_mm = dist_mm - Self::MAX_DISTANCE_MM;
                let max_dist_mm = dist_mm + Self::MAX_DISTANCE_MM;
                let start = bisect_left_by(&self.point_distances, |p| p.0.cmp(&min_dist_mm));
                for it in self.point_distances[start..].iter() {
                    if it.0 > max_dist_mm {
                        break;
                    }
                    let p1 = &self.points[it.1 as usize];
                    let p2 = &self.points[it.2 as usize];
                    if let Some(tr) = make_transform(s1, p1, s2, p2) {
                        let score = self.score_transform(&tr);
                        if score.better_than(&best_score) {
                            best_transform = tr;
                            best_score = score;
                            // println!("best_score={:?}", score)
                        }
                    }
                    if let Some(tr) = make_transform(s1, p2, s2, p1) {
                        let score = self.score_transform(&tr);
                        if score.better_than(&best_score) {
                            best_transform = tr;
                            best_score = score;
                            // println!("best_score={:?}", score);
                        }
                    }
                }
            }
        }

        // At least half of symbols should match, otherwise we have
        // no trust in the result.
        if best_score.num_matches >= num_symbols / 2 {
            let result = self.refine_transform(&best_transform);
            let score = self.score_transform(&result);
            Some((result, score))
        } else {
            None
        }
    }

    // TODO: Refactor the common logic of score_transform() and refine_transform()
    // into a single place. Maybe an iterator over (&Symbol, &Point) tuples?
    fn score_transform(&self, gt: &GeoTransform) -> Score {
        // Must be kept in sync with refine_transform().
        let max_dist_m = (Self::MAX_DISTANCE_MM as f64) / 1000.0;
        let max_dist_sq = max_dist_m * max_dist_m;

        let mut num_matches = 0;
        let mut sum_dist = 0.0;
        for sym in self.symbols.iter() {
            let (x, y) = sym.project(gt);
            for (_p, dist_sq) in self
                .point_tree
                .nearest_neighbor_iter_with_distance_2(&[x, y])
            {
                if dist_sq > max_dist_sq {
                    break;
                }

                // TODO: Once we have symbol types, guard the
                // following by "if sym.type == p.type {...}".
                if true {
                    num_matches += 1;
                    // In benchmarks, computing the square root did
                    // not lead to any measureable difference in speed.
                    // Since the average distance of matches is
                    // a meaningful metric, whereas the summed squares
                    // are more difficult to interpret over a dataset
                    // with wildly varying num_matches, we do compute
                    // the square root here.
                    sum_dist += dist_sq.sqrt();
                    break;
                }
            }
        }

        // We always have at least two matches due to the way how
        // the candidate transforms are constructed, so the division
        // by num_matches is safe.
        Score {
            num_matches,
            avg_dist: sum_dist / (num_matches as f64),
        }
    }

    fn refine_transform(&self, gt: &GeoTransform) -> GeoTransform {
        // Must be kept in sync with score_transform().
        let max_dist_m = (Self::MAX_DISTANCE_MM as f64) / 1000.0;
        let max_dist_sq = max_dist_m * max_dist_m;

        let mut gcps = Vec::with_capacity(self.symbols.len());
        for sym in self.symbols.iter() {
            let (x, y) = sym.project(gt);
            for (p, dist_sq) in self
                .point_tree
                .nearest_neighbor_iter_with_distance_2(&[x, y])
            {
                if dist_sq > max_dist_sq {
                    break;
                }

                // TODO: Once we have symbol types, guard the
                // following by "if sym.type == p.type {...}".
                if true {
                    gcps.push(make_gcp(sym, &self.points[p.data]));
                    break;
                }
            }
        }

        let mut result = [0.0; 6];
        let ok = unsafe { GDALGCPsToGeoTransform(2, gcps.as_ptr(), result.as_mut_ptr(), 0) != 0 };
        if !ok {
            return *gt;
        }

        // When GDAL computes a transform from *all* matching Ground Control Points,
        // the result can be worse than our own, initial estimate that was just
        // computed from two points. We threfore check the quality of the refinement,
        // and return the refined, GDAL-estimated transform only if it s actually
        // better than our own.
        let score = self.score_transform(&result);
        let unrefined_score = self.score_transform(gt);
        if score.better_than(&unrefined_score) {
            result
        } else {
            *gt
        }
    }
}

impl Symbol {
    pub fn distance_mm(&self, other: &Self, meters_per_pixel: f64) -> i32 {
        let dx = (self.x - other.x) * meters_per_pixel;
        let dy = (self.y - other.y) * meters_per_pixel;
        let dist_sq = dx * dx + dy * dy;

        (dist_sq.sqrt() * 1000.0 + 0.5) as i32
    }

    pub fn project(&self, gt: &GeoTransform) -> (f64, f64) {
        // https://gdal.org/en/latest/tutorials/geotransforms_tut.html
        // X_geo = GT(0) + X_pixel * GT(1) + Y_line * GT(2)
        // Y_geo = GT(3) + X_pixel * GT(4) + Y_line * GT(5)
        let x = gt[0] + self.x * gt[1] + self.y * gt[2];
        let y = gt[3] + self.x * gt[4] + self.y * gt[5];
        (x, y)
    }
}

// Calls GDAL to open a TIFF file, returning a GDAL Dataset.
fn tiff_open(img: &PathBuf, page: u32) -> Result<Dataset, Box<dyn Error>> {
    // GDAL supports opening multi-page TIFF, but the API is a little arcane:
    // instead of making this an explicit API parameter, the desired page
    // needs to be passed as part of the file name, using a special GTIFF_DIR
    // prefix.
    let path = <PathBuf as Clone>::clone(img);
    let path = path.into_os_string().into_string();
    if path.is_err() {
        return Err("cannot convert OSString to Unicode string".into());
    };
    let path = path.unwrap();
    let path = format!("GTIFF_DIR:{page}:{path}");
    let img = Dataset::open(PathBuf::from(path))?;
    Ok(img)
}

// Returns the resolution of the passed image in pixels per inch.
pub fn image_resolution_dpi(img: &PathBuf, page: u32) -> Result<f64, Box<dyn Error>> {
    let img = tiff_open(img, page)?;
    const DEFAULT_DOMAIN: &str = ""; // namespace for GDAL metadata
    let xres = img.metadata_item("TIFFTAG_XRESOLUTION", DEFAULT_DOMAIN);
    let yres = img.metadata_item("TIFFTAG_YRESOLUTION", DEFAULT_DOMAIN);
    let unit = img.metadata_item("TIFFTAG_RESOLUTIONUNIT", DEFAULT_DOMAIN);
    img.close()?;

    let unit_scale = match unit.unwrap_or(String::from("")).as_str() {
        "2 (pixels/inch)" => 1.0,
        "3 (pixels/cm)" => 2.54,
        _ => Err("image does not specify a resolution unit")?,
    };

    let xres = xres.unwrap_or(String::from("72")).parse::<f64>()?;
    let yres = yres.unwrap_or(String::from("72")).parse::<f64>()?;
    if xres != yres {
        return Err("x and y resolution must be the same".into());
    }

    Ok(xres * unit_scale)
}

pub fn write_geotiff(
    img: PathBuf,
    page: u32,
    gt: &GeoTransform,
    out: PathBuf,
) -> Result<(), Box<dyn Error>> {
    // "cog" = Cloud-Optimized GeoTIFF
    // https://gdal.org/en/latest/drivers/raster/cog.html
    let driver = DriverManager::get_driver_by_name("cog")?;
    let mut opts = RasterCreationOptions::new();
    opts.set_name_value("NUM_THREADS", "ALL_CPUS")?;
    opts.set_name_value("STATISTICS", "YES")?;
    opts.set_name_value("PREDICTOR", "YES")?;
    opts.set_name_value("COMPRESS", "DEFLATE")?;

    let img = tiff_open(&img, page)?;
    let mut out = img.create_copy(&driver, out, &opts)?;
    out.set_projection("epsg:2056")?; // epsg.io/2056 = Swiss CH1903+/LV95
    out.set_geo_transform(gt)?;

    out.close()?;
    img.close()?;
    Ok(())
}

fn read_points(path: PathBuf) -> Result<Vec<Point>, Box<dyn Error>> {
    let mut points = Vec::new();
    let mut reader = Reader::from_reader(File::open(path)?);
    for rec in reader.deserialize() {
        let p: Point = rec?;
        points.push(p);
    }
    sort_points(&mut points);
    Ok(points)
}

fn sort_points(points: &mut [Point]) {
    if points.len() < 2 {
        return;
    }

    let x0 = points
        .iter()
        .min_by(|a, b| a.x.partial_cmp(&b.x).unwrap())
        .unwrap()
        .x;
    let y0 = points
        .iter()
        .min_by(|a, b| a.y.partial_cmp(&b.y).unwrap())
        .unwrap()
        .y;
    points.sort_by(|a, b| {
        let da = (a.x - x0) * (a.x - x0) + (a.y - y0) * (a.y - y0);
        let db = (b.x - x0) * (b.x - x0) + (b.y - y0) * (b.y - y0);
        da.partial_cmp(&db).unwrap()
    });
}

fn read_symbols(path: PathBuf) -> Result<Vec<Symbol>, Box<dyn Error>> {
    let mut symbols = Vec::new();
    let mut reader = Reader::from_reader(File::open(path)?);
    for rec in reader.deserialize() {
        let p: Symbol = rec?;
        symbols.push(p);
    }
    sort_symbols(&mut symbols);
    Ok(symbols)
}

fn sort_symbols(symbols: &mut [Symbol]) {
    if symbols.len() < 2 {
        return;
    }

    let x0 = symbols
        .iter()
        .min_by(|a, b| a.x.partial_cmp(&b.x).unwrap())
        .unwrap()
        .x;
    let y0 = symbols
        .iter()
        .min_by(|a, b| a.y.partial_cmp(&b.y).unwrap())
        .unwrap()
        .y;
    symbols.sort_by(|a, b| {
        let da = (a.x - x0) * (a.x - x0) + (a.y - y0) * (a.y - y0);
        let db = (b.x - x0) * (b.x - x0) + (b.y - y0) * (b.y - y0);
        da.partial_cmp(&db).unwrap()
    });
}

fn make_gcp(s: &Symbol, p: &Point) -> GDAL_GCP {
    GDAL_GCP {
        pszId: std::ptr::null_mut(),
        pszInfo: std::ptr::null_mut(),
        dfGCPPixel: s.x,
        dfGCPLine: s.y,
        dfGCPX: p.x,
        dfGCPY: p.y,
        dfGCPZ: 0.0,
    }
}

// Create a GDAL transform. Returns None for some corner cases,
// such as two symbols having the exact same y coordinate.
fn make_transform(s1: &Symbol, p1: &Point, s2: &Symbol, p2: &Point) -> Option<GeoTransform> {
    let mut result = [0.0; 6];
    let gcps = [make_gcp(s1, p1), make_gcp(s2, p2)];
    let ok = unsafe { GDALGCPsToGeoTransform(2, gcps.as_ptr(), result.as_mut_ptr(), 0) != 0 };
    if ok {
        Some(result)
    } else {
        None
    }
}

fn build_point_distances(pts: &[Point]) -> PointDistances {
    let cap = (pts.len() * (pts.len() - 1)) / 2;
    let mut result = PointDistances::with_capacity(cap);
    for i in 0..pts.len() {
        for j in i + 1..pts.len() {
            let dx = pts[i].x - pts[j].x;
            let dy = pts[i].y - pts[j].y;
            let dist_sq = dx * dx + dy * dy;
            let dist = (dist_sq.sqrt() * 1000.0 + 0.5) as i32;
            result.push((dist, i as u32, j as u32));
        }
    }
    result.sort_by_key(|a| a.0);
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn can_sort_symbols() {
        let mut syms = vec![Symbol { x: 4., y: 3. }, Symbol { x: 1., y: 2. }];
        sort_symbols(&mut syms);
        let x_iter = syms.iter().map(|p| format!("{:.0}", p.x));
        let xs = x_iter.collect::<Vec<String>>().join(",");
        assert_eq!(xs, "1,4");
    }

    #[test]
    fn can_sort_points() {
        let mut points = vec![
            Point {
                id: "A".to_string(),
                x: 4.,
                y: 2.,
            },
            Point {
                id: "B".to_string(),
                x: 77.,
                y: 99.0,
            },
            Point {
                id: "C".to_string(),
                x: 11.,
                y: 33.0,
            },
        ];
        sort_points(&mut points);
        let id_iter = points.iter().map(|p| p.id.as_str());
        let ids = id_iter.collect::<Vec<&str>>().join(",");
        assert_eq!(ids, "A,C,B");
    }

    #[test]
    fn can_make_gcp() {
        let gcp = make_gcp(
            &Symbol {
                x: 1723.25,
                y: 3282.0,
            },
            &Point {
                id: "HGF4845".to_string(),
                x: 2679878.122,
                y: 1251033.022,
            },
        );

        // assert_eq!(gcp.id, "HGF4845");
        // assert_eq!(gcp.info, "");
        assert_eq!(gcp.dfGCPPixel, 1723.25);
        assert_eq!(gcp.dfGCPLine, 3282.0);
        assert_eq!(gcp.dfGCPX, 2679878.122);
        assert_eq!(gcp.dfGCPY, 1251033.022);
        assert_eq!(gcp.dfGCPZ, 0.0);
    }

    #[test]
    fn can_build_point_distances() {
        let points = vec![
            Point {
                id: "HGF4845".to_string(),
                x: 2679878.122,
                y: 1251033.022,
            },
            Point {
                id: "HGC4303".to_string(),
                x: 2678592.752,
                y: 1251921.679,
            },
            Point {
                id: "HGB4303".to_string(),
                x: 2678599.342,
                y: 1251903.290,
            },
        ];
        let dist = build_point_distances(&points);
        assert_eq!(dist.len(), 3);
        assert_eq!(dist[0], (19534, 1, 2));
        assert_eq!(dist[1], (1546818, 0, 2));
        assert_eq!(dist[2], (1562654, 0, 1));
    }

    #[test]
    fn can_find_image_resolution_dpi() {
        let mut path = std::env::current_dir().unwrap();
        path.push("testdata");

        path.push("HG3099.tif");
        assert_eq!(image_resolution_dpi(&path, 1).unwrap(), 300.0);

        path.set_file_name("500_pixels_per_cm.tif");
        assert_eq!(image_resolution_dpi(&path, 1).unwrap(), 1270.0);

        path.set_file_name("different_x_and_y_resolution.tif");
        assert_eq!(image_resolution_dpi(&path, 1).is_err(), true);
    }
}
