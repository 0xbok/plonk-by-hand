mod elliptic_curve;

use crate::elliptic_curve::{ECurve, Point};

fn main() {
    let ec = ECurve{ field: 101 };
    let p = Point{ x: 1, y: 2};
    let mut q = ec.mul(&p, 2);
    println!("{:?}", q);

    for _ in 3..17 {
        q = ec.add(&q, &p);
        println!("{:?}", q);
    }

    // println!("{:?}",ec.mul(&p, 2));
}
