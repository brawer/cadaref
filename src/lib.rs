use csv::Reader;
use gdal::Gcp;
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

impl Matcher {
    pub fn new(points: PathBuf, symbols: PathBuf) -> Result<Self, Box<dyn Error>> {
        let points = read_points(points)?;
        let mut point_tree = RTree::new();
        for (i, p) in points.iter().enumerate() {
            point_tree.insert(PointNode::new([p.x, p.y], i));
        }

        let point_distances = build_point_distances(&points);
        // TODO: Build list of distances between points;
        // sort by distance; build look-up table from
        // int(distance) -> index in point-distance table;
        // write fn points_in_distance(dist: f32) -> iter(pointA, pointB, dist).

        let symbols = read_symbols(symbols)?;

        let matcher = Matcher {
            points,
            point_distances,
            point_tree,
            symbols,
        };
        Ok(matcher)
    }

    pub fn find_matches(&self, meters_per_pixel: f64) {
        for i in 0..self.symbols.len() {
            for j in i + 1..self.symbols.len() {
                let p = &self.symbols[i];
                let q = &self.symbols[j];
                let dist_mm = p.distance_mm(q, meters_per_pixel);
                //let dx = (p.x - q.x) * meters_per_pixel;
                //let dy = (p.y - q.y) * meters_per_pixel;
                //let dist_sq = dx*dx + dy*dy;
                //let dist_mm = (dist_sq.sqrt() * 1000.0 + 0.5) as i32;

                println!(
                    "i={:?} j={:?} dist_mm={dist_mm}",
                    self.symbols[i], self.symbols[j]
                );
            }
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

fn make_gcp(s: &Symbol, p: &Point) -> Gcp {
    Gcp {
        id: p.id.clone(),
        info: String::new(),
        pixel: s.x,
        line: s.y,
        x: p.x,
        y: p.y,
        z: 0.,
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
        assert_eq!(gcp.id, "HGF4845");
        assert_eq!(gcp.info, "");
        assert_eq!(gcp.pixel, 1723.25);
        assert_eq!(gcp.line, 3282.0);
        assert_eq!(gcp.x, 2679878.122);
        assert_eq!(gcp.y, 1251033.022);
        assert_eq!(gcp.z, 0.0);
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
