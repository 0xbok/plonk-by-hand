mod elliptic_curve;

use crate::elliptic_curve::{ECurve, Point};

fn main() {
    let ec = ECurve{ field: 101 };
    let p = Point{ x: 1, y: 2};

    println!("{:?}", ec.neg(&ec.double(&p)));
}
