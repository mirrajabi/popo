use std::time::Instant;

use popo::poisson_disc_sampling::sample;
use popo::vectors::Vec2;

fn main() {
    let polygons = seed_polygons();
    loop {
        for polygon in polygons.clone() {
            let polygon: Vec<Vec2> = polygon.iter().map(|v| Vec2::new(v.x, v.y)).collect();
            let start = Instant::now();
            let samples: Vec<Vec2> = sample(polygon.clone(), 0.1, 30, Some(0.), None).collect();
            let elapsed = start.elapsed().as_millis();
            println!("Took: {elapsed}ms to generate {} samples", samples.len());
        }
    }
}

fn seed_polygons() -> Vec<Vec<Vec2>> {
    vec![vec![
        Vec2::new(0., 0.),
        Vec2::new(0., -20.),
        Vec2::new(2., -20.),
        Vec2::new(2., -30.),
        Vec2::new(18., -28.),
        Vec2::new(22., -22.),
        Vec2::new(26., -25.),
        Vec2::new(28., -15.),
        Vec2::new(25., -10.),
        Vec2::new(20., -8.),
        Vec2::new(18., -4.),
        Vec2::new(20., 0.),
        Vec2::new(10., 0.),
        Vec2::new(10., 6.),
        Vec2::new(12., 6.),
        Vec2::new(12., 16.),
        Vec2::new(8., 18.),
        Vec2::new(10., 26.),
        Vec2::new(0., 28.),
        Vec2::new(-4., 24.),
        Vec2::new(-10., 30.),
        Vec2::new(-16., 26.),
        Vec2::new(-14., 20.),
        Vec2::new(-22., 18.),
        Vec2::new(-24., 12.),
        Vec2::new(-28., 8.),
        Vec2::new(-26., 4.),
        Vec2::new(-20., 2.),
        Vec2::new(-20., -6.),
        Vec2::new(-14., -6.),
        Vec2::new(-10., -2.),
        Vec2::new(-10., 0.),
        Vec2::new(0., 0.),
    ]]
}
