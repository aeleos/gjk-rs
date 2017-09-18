#![allow(dead_code)]

#[cfg(test)]
mod tests {

    extern crate rand;

    use super::*;

    pub const EPSILON: f64 = 2.2204460492503131e-16f64;

    fn perturbation() -> f64 {
        (rand::random::<u32>() as f64 / u32::max_value() as f64) * EPSILON * 100.0 *
            if rand::random::<u8>() % 2 == 1 {
                1.0
            } else {
                -1.0
            }
    }

    fn jostle(vec: &mut Vec2) -> &Vec2 {
        vec.x = vec.x + perturbation();
        vec.y = vec.y + perturbation();
        vec
    }

    fn jostle_copy(vec: &Vec2) -> Vec2 {
        Vec2::new(vec.x + perturbation(), vec.y + perturbation())
    }




    #[test]
    fn test_random_detections() {
        let vertices1 = [
            Vec2::new(4.0, 11.0),
            Vec2::new(5.0, 5.0),
            Vec2::new(9.0, 9.0),
        ];
        let vertices2 = [
            Vec2::new(4.0, 11.0),
            Vec2::new(5.0, 5.0),
            Vec2::new(9.0, 9.0),
        ];

        let mut num_no_collisions = 0;

        for iter in 0..100 {
            let a = vertices1
                .iter()
                .map(|&v| jostle_copy(&v))
                .collect::<Vec<Vec2>>();

            let b = vertices2
                .iter()
                .map(|&v| jostle_copy(&v))
                .collect::<Vec<Vec2>>();

            if gjk(&a[..], &b[..]) == 1 {
                num_no_collisions += 1;
            }
        }

        assert!(num_no_collisions > 0);
    }
}


#[derive(Debug, Copy, Clone)]
struct Vec2 {
    x: f64,
    y: f64,
}

impl Vec2 {
    fn new(x: f64, y: f64) -> Vec2 {
        Vec2 { x: x, y: y }
    }

    fn subtract(&mut self, a: &Vec2) -> &mut Vec2 {
        self.x -= a.x;
        self.y -= a.y;
        self
    }

    fn subtract_copy(a: &Vec2, b: &Vec2) -> Vec2 {
        Vec2::new(a.x - b.x, a.y - b.y)
    }

    fn negate(&mut self) -> &mut Vec2 {
        self.x = -self.x;
        self.y = -self.y;
        self
    }

    fn negate_copy(&self) -> Vec2 {
        Vec2::new(-self.x, -self.y)
    }

    fn perpendicular(&mut self) -> &mut Vec2 {
        let temp_x = self.x;
        self.x = self.y;
        self.y = -temp_x;
        self
    }

    fn perpendicular_copy(&self) -> Vec2 {
        Vec2::new(self.y, -self.x)
    }

    fn dot_product(&self, other_vec: &Vec2) -> f64 {
        (self.x * other_vec.x) + (self.y * other_vec.y)
    }

    fn length_squared(&self) -> f64 {
        (self.x * self.x) + (self.y * self.y)
    }

    /// Triple product expansion is used to calculate perpendicular normal vectors
    /// which kinda 'prefer' pointing towards the Origin in Minkowski space
    fn triple_product(a: &Vec2, b: &Vec2, c: &Vec2) -> Vec2 {
        let ac = a.dot_product(c);
        let bc = b.dot_product(c);
        Vec2::new(b.x * ac - a.x * bc, b.y * ac - a.y * bc)
    }

    fn average_point(vertices: &[Vec2]) -> Vec2 {
        assert!(vertices.len() > 0);
        let mut avg = Vec2::new(0.0 as f64, 0.0 as f64);
        for vert in vertices.iter() {
            avg.x += vert.x;
            avg.y += vert.y;
        }

        avg.x /= vertices.len() as f64;
        avg.y /= vertices.len() as f64;
        avg
    }

    fn index_of_furthest_point(&self, vertices: &[Vec2]) -> usize {
        assert!(vertices.len() > 0);
        let mut max_product = self.dot_product(&vertices[0]);
        let mut idx: usize = 0;
        for (i, vert) in vertices.iter().enumerate() {
            let product = self.dot_product(vert);
            if product > max_product {
                max_product = product;
                idx = i;
            }
        }
        idx
    }

    fn minkowski_sum(&mut self, vertices1: &[Vec2], vertices2: &[Vec2]) -> Vec2 {
        let i = self.index_of_furthest_point(vertices1);
        let j = self.negate_copy().index_of_furthest_point(vertices2);
        Vec2::subtract_copy(&vertices1[i], &vertices2[j])
    }
}


fn gjk(vertices1: &[Vec2], vertices2: &[Vec2]) -> isize {
    let mut index: usize = 0;
    // let a, b, c, d, ao, ab, ac, abperp, acperp, simplex[3]: Vec2;

    let position1 = Vec2::average_point(vertices1);
    let position2 = Vec2::average_point(vertices2);

    let mut d = Vec2::subtract_copy(&position1, &position2);

    if (d.x == 0.0) && (d.y == 0.0) {
        d.x = 1.0;
    }

    let mut simplex = [Vec2 { x: 0.0, y: 0.0 }; 3];

    let mut a = d.minkowski_sum(vertices1, vertices2);

    simplex[0] = a;

    if a.dot_product(&d) <= 0.0 {
        return 0;
    }

    d = a.negate_copy();

    let (mut b, mut c, mut ao, mut ab, mut ac, mut abperp, mut acperp);

    loop {
        index += 1;

        simplex[index] = d.minkowski_sum(vertices1, vertices2);

        a = simplex[index];

        if a.dot_product(&d) <= 0.0 {
            return 0;
        }

        ao = a.negate_copy();

        if index < 2 {
            b = simplex[0];
            ab = Vec2::subtract_copy(&b, &a);
            d = Vec2::triple_product(&ab, &ao, &ab);
            if d.length_squared() == 0.0 {
                d = ab.perpendicular_copy();
            }
            continue;
        }

        b = simplex[1];
        c = simplex[0];
        ab = Vec2::subtract_copy(&b, &a);
        ac = Vec2::subtract_copy(&c, &a);

        acperp = Vec2::triple_product(&ab, &ac, &ac);

        if acperp.dot_product(&ao) >= 0.0 {
            d = acperp;
        } else {
            abperp = Vec2::triple_product(&ac, &ab, &ab);

            if abperp.dot_product(&ao) < 0.0 {
                return 1;
            }

            simplex[0] = simplex[1];

            d = abperp;
        }

        simplex[1] = simplex[2];
        index -= 1;
    }
}
