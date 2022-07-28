// elliptic curve y^2 = x^3 + ax + b where a = 0

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

        return Self::gcd (a, b-a);
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
        return r;
    }

    // p -> 2p
    pub fn double(&self, p: &Point) -> Point {
        let m_a = 3 * p.x * p.x;
        let m_b = 2 * p.y;
        let g = Self::gcd(m_a, m_b);
        let m_a = m_a / g;
        let m_b = m_b / g;

        // https://math.stackexchange.com/a/586613
        // https://cp-algorithms.com/algebra/phi-function.html
        let m_b_mul_inv = self.pow(m_b, self.field - 2);
        let m = self.modular((m_a * m_b_mul_inv) as i32);

        return Point {
            x: self.modular((m*m) as i32 - (2*p.x) as i32),
            y: self.modular(m as i32 * (3*p.x as i32 - (m as i32 * m as i32)) - p.y as i32)
        }
    }
}