use bisection::bisect_left_by;
use csv::Reader;
use gdal::GeoTransform;
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
struct Score {
    num_matches: usize,
    sum_dist_sq: f64,
}

impl Score {
    fn better_than(&self, other: &Self) -> bool {
        if self.num_matches > other.num_matches {
            return true;
        }
        if self.num_matches < other.num_matches {
            return false;
        }
        self.sum_dist_sq < other.sum_dist_sq
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

    pub fn find_matches(&self, meters_per_pixel: f64) {
        let num_symbols = self.symbols.len();
        println!("start");

        let mut best_transform = [0.0; 6];
        let mut best_score = Score {
            num_matches: 0,
            sum_dist_sq: 0.0,
        };

        const MAX_MATCHES: usize = 10;
        for i in 0..num_symbols.max(MAX_MATCHES) {
            let mut start = i + 1;
            if num_symbols - start > MAX_MATCHES {
                start = num_symbols - MAX_MATCHES;
            }

            for j in start..num_symbols {
                let s1 = &self.symbols[i];
                let s2 = &self.symbols[j];
                let dist_mm = s1.distance_mm(s2, meters_per_pixel);
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

        let refined = self.refine_transform(&best_transform);
        let refined_score = self.score_transform(&refined);
        println!("best_score={:?} for {:?}", best_score, best_transform);
        println!("final_score={:?} for {:?}", refined_score, refined);
    }

    // TODO: Refactor the common logic of score_transform() and refine_transform()
    // into a single place. Maybe an iterator over (&Symbol, &Point) tuples?
    fn score_transform(&self, gt: &GeoTransform) -> Score {
        // Must be kept in sync with refine_transform().
        let max_dist_m = (Self::MAX_DISTANCE_MM as f64) / 1000.0;
        let max_dist_sq = max_dist_m * max_dist_m;

        let mut num_matches = 0;
        let mut sum_dist_sq = 0.0;
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
                    sum_dist_sq += dist_sq;
                    break;
                }
            }
        }

        Score {
            num_matches,
            sum_dist_sq,
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
        if ok {
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
}
