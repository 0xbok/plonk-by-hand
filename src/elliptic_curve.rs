// elliptic curve y^2 = x^3 + ax + b where a = 0
// b is automatically determined from a point on the curve

pub struct ECurve {
    pub field: u32,
}

#[derive(Debug, Clone, Copy)]
#[derive(PartialEq, Eq)]
pub struct Point {
    pub x: u32,
    pub y: u32,
}

// each function assumes the point's coordinate is in [0, field-1]
impl ECurve {
    // p -> -p
    pub fn neg(&self, p: &Point) -> Point {
        Point {
            x: p.x,
            y: self.field - p.y
        }
    }

    fn gcd(a: u32, b: u32) -> u32 {
        if a == 0 {
            return b;
        } else if a > b {
            return Self::gcd(b, a);
        }

        Self::gcd(a, b-a)
    }

    fn modular(&self, a: i32) -> u32 {
        if a >= 0 {
            return a as u32 % self.field;
        } else {
            return self.field - (-a as u32 % self.field);
        }
    }

    // b^e mod self.field
    // https://en.wikipedia.org/wiki/Modular_exponentiation#Implementation_in_Lua
    fn pow(&self, mut b: u32, mut e: u32) -> u32 {
        if self.field == 1 {
            return 0;
        }
        let mut r = 1;
        while e > 0 {
            if e & 1 == 1 {
                r = (r * b) % self.field;
            }
            b = (b*b) % self.field;
            e >>= 1;
        }
        r
    }

    // convert fraction to a field element
    fn frac_to_element(&self, a: u32, b: u32) -> u32 {
        let g: u32 = Self::gcd(a, b);
        let a = a / g;
        let b = b / g;

        // https://math.stackexchange.com/a/586613
        // https://cp-algorithms.com/algebra/phi-function.html
        let b_mul_inv = self.pow(b, self.field - 2);
        self.modular((a * b_mul_inv) as i32)
    }

    // p -> 2p
    fn double(&self, p: &Point) -> Point {
        // m = m_a/m_b
        let m_a = 3 * p.x * p.x;
        let m_b = 2 * p.y;
        let m = self.frac_to_element(m_a, m_b);

        Point {
            x: self.modular((m*m) as i32 - (2*p.x) as i32),
            y: self.modular(m as i32 * (3*p.x as i32 - (m as i32 * m as i32)) - p.y as i32)
        }
    }

    pub fn add(&self, p: &Point, q: &Point) -> Point {
        assert_ne!(p.x, q.x);

        let lam_a = self.modular(q.y as i32 - p.y as i32);
        let lam_b = self.modular(q.x as i32 - p.x as i32);

        let lam = self.frac_to_element(lam_a, lam_b);

        let r_x = self.modular((lam*lam) as i32 - (p.x + q.x) as i32);

        Point {
            x: r_x,
            y: self.modular((lam*p.x) as i32 - (lam*r_x + p.y) as i32),
        }
    }

    // p -> n*p
    pub fn mul(&self, p: &Point, n: u32) -> Point {
        assert_eq!(n & 1, 0);
        if n == 2 {
            return self.double(p);
        }
        let r = self.mul(p, n/2);
        self.double(&r)
    }
}