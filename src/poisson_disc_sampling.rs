use crate::{
    polygons::{Polygon, offset},
    vectors::Vec2,
};
use core::panic;
use std::{f32, iter};

/// A function that places discs in a polygon given the ordered polygon vertices and the radius of each disc,
/// in a way that no two disc will have a collision. It does not care about the winding order of given vertices,
/// so they can be passed in either CW or CCW order.
///
/// # Arguments
/// - [`polygon`] A vector of `Vec2` points representing the vertices of the polygon. It can be presented in either CW or CCW order.
/// - [`r`]: The radius of the discs to be placed.
/// - [`max_attempts`]: The maximum number of attempts to place a disc before giving up. Typically 30 is advised.
/// - [`padding`]: How much padding should be considered along the edges. A positive padding value indicates inwards
///     offset, and a negative padding indicates outwards offset.
/// - [`start_point`]: An optional starting point for the first disc. If it is not inside the polygon or is not provided,
///     a random point within the polygon will be chosen.
///
/// # Returns
///
/// A lazy iterator over `Vec2` points sampled within the polygon. Each point is guaranteed
/// to be at least `r` distance from every other point.
///
/// # Examples
///
/// ```rust
/// use popo::poisson_disc_sampling::sample;
/// use popo::vectors::Vec2;
///
/// let polygon = vec![
///     Vec2::new(0.0, 0.0),
///     Vec2::new(1.0, 0.0),
///     Vec2::new(1.0, 1.0),
///     Vec2::new(0.0, 1.0),
/// ];
/// let samples: Vec<Vec2> = sample(polygon, 0.1, 30, Some(0.01), None).collect();
/// assert!(!samples.is_empty());
/// println!("Generated {} samples", samples.len());
/// for point in samples {
///     println!("{:?}", point);
/// }
/// ```
///
/// # The Poisson Disc algorithm
///
/// The algorithm takes as input the extent of the sample domain in Rn,
/// the minimum distance r between samples, and a constant k as the limit
/// of samples to choose before rejection in the algorithm (typically k = 30).
///
/// Step 0. Initialize an n-dimensional background grid for storing
/// samples and accelerating spatial searches. We pick the cell size to
/// be bounded by r/√n, so that each grid cell will contain at most
/// one sample, and thus the grid can be implemented as a simple n dimensional
/// array of integers:
/// - the default −1 indicates no sample,
/// - a non-negative integer gives the index of the sample located in a cell.
///
/// Step 1. Select the initial sample, x0, randomly chosen uniformly
/// from the domain. Insert it into the background grid, and initialize
/// the “active list” (an array of sample indices) with this index (zero).
///
/// Step 2. While the active list is not empty, choose a random index
/// from it (say i). Generate up to k points chosen uniformly from the
/// spherical annulus between radius r and 2r around xi. For each
/// point in turn, check if it is within distance r of existing samples
/// (using the background grid to only test nearby samples). If a point
/// is adequately far from existing samples, emit it as the next sample
/// and add it to the active list. If after k attempts no such point is
/// found, instead remove i from the active list.
///
/// **Source:** https://www.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf
///
/// # Notes
/// - The current implementation does not support polygons with holes.
/// - The current implementation is single-threaded
pub fn sample(
    mut polygon: Polygon,
    r: f32,
    max_attempts: u16,
    padding: Option<f32>,
    start_point: Option<Vec2>,
) -> Box<dyn Iterator<Item = Vec2>> {
    if polygon.len() < 3 {
        return Box::new(iter::empty());
    }

    let first = polygon.first().unwrap();
    let last = polygon.last().unwrap();
    if first.x != last.x || first.y != last.y {
        polygon.push(polygon[0]);
    }

    if let Some(p) = padding {
        // Offset function has an inverted logic for positive and negative padding, so we need to invert the padding value.
        polygon = offset(polygon, -p);
    }

    let n_dimensions = 2f32;
    let cell_size = r / n_dimensions.sqrt();
    let mut polygon_clone = polygon.clone();
    polygon_clone.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap());
    let min_x = polygon_clone.first().unwrap().x;
    let max_x = polygon_clone.last().unwrap().x;
    polygon_clone = polygon.clone();
    polygon_clone.sort_by(|a, b| a.y.partial_cmp(&b.y).unwrap());
    let min_y = polygon_clone.first().unwrap().y;
    let max_y = polygon_clone.last().unwrap().y;

    let w = (max_x - min_x).abs();
    let h = (max_y - min_y).abs();
    let sample_w = (w / cell_size).ceil() as usize;
    let sample_h = (h / cell_size).ceil() as usize;
    let size = (sample_w as usize).checked_mul(sample_h as usize).expect(
        "Buffer size overflow. The result of (sampleW * sampleH) is too big to fit in usize.",
    );
    let mut grid = vec![0; size];

    let mut points = Vec::<Vec2>::new();
    let mut active_list = Vec::<Vec2>::new();

    let mut rng = fastrand::Rng::new();

    let mut initial_point = if let Some(sp) = start_point {
        sp
    } else {
        Vec2 {
            x: min_x + rng.f32() * (max_x - min_x),
            y: min_y + rng.f32() * (max_y - min_y),
        }
    };
    while !is_in_polygon(initial_point, &polygon) {
        initial_point = Vec2 {
            x: (min_x + rng.f32() * (max_x - min_x)),
            y: (min_y + rng.f32() * (max_y - min_y)),
        };
    }

    active_list.push(initial_point);

    let iterator = iter::from_fn(move || {
        let index_fn = |i: usize, j: usize| -> usize { j * sample_w + i };

        loop {
            if active_list.len() == 0 {
                return None;
            }

            let active_point_index = rng.usize(0..active_list.len());
            let active_point = active_list[active_point_index];

            let mut attempt = 0;

            while attempt < max_attempts {
                let angle = rng.f32() * f32::consts::PI * 2.;
                let radius = r + rng.f32() * (2.0 - r) * r;
                let direction = Vec2::from_angle(angle);
                let candidate = active_point + direction * radius;

                let cell_x = ((candidate.x - min_x) / cell_size).floor() as usize;
                let cell_y = ((candidate.y - min_y) / cell_size).floor() as usize;
                if is_valid(
                    candidate,
                    &polygon,
                    min_x,
                    max_x,
                    min_y,
                    max_y,
                    sample_w,
                    sample_h,
                    cell_x,
                    cell_y,
                    r,
                    &points,
                    &grid,
                    index_fn,
                ) {
                    points.push(candidate);
                    active_list.push(candidate);
                    let index = index_fn(cell_x, cell_y);
                    grid[index] = points.len();
                    return Some(candidate);
                }

                attempt += 1;
            }

            active_list.remove(active_point_index);

            if active_list.len() > size {
                panic!(
                    "Algorithm fault! The algorithm has produced more points than there could possibly be. Something is terribly wrong!"
                )
            }
        }
    });
    Box::new(iterator)
}

fn is_in_polygon(candidate: Vec2, polygon: &[Vec2]) -> bool {
    let n_vert = polygon.len();
    let mut j = n_vert - 1;
    let mut inside = false;
    for i in 0..n_vert {
        let vert_i = polygon[i];
        let vert_j = polygon[j];
        if ((vert_i.y > candidate.y) != (vert_j.y > candidate.y))
            && (candidate.x
                < (vert_j.x - vert_i.x) * (candidate.y - vert_i.y) / (vert_j.y - vert_i.y)
                    + vert_i.x)
        {
            inside = !inside;
        }
        j = i;
    }

    return inside;
}

fn is_valid<IndexFn>(
    candidate: Vec2,
    polygon: &[Vec2],
    min_x: f32,
    max_x: f32,
    min_y: f32,
    max_y: f32,
    sample_w: usize,
    sample_h: usize,
    cell_x: usize,
    cell_y: usize,
    r: f32,
    points: &[Vec2],
    grid: &[usize],
    index_fn: IndexFn,
) -> bool
where
    IndexFn: Fn(usize, usize) -> usize,
{
    if candidate.x < min_x || candidate.x >= max_x || candidate.y < min_y || candidate.y >= max_y {
        return false;
    }

    let range_x_start = cell_x.checked_sub(2).unwrap_or(0);
    let range_x_end = (cell_x + 2).min(sample_w);
    let range_y_start = cell_y.checked_sub(2).unwrap_or(0);
    let range_y_end = (cell_y + 2).min(sample_h);
    for y in range_y_start..range_y_end {
        for x in range_x_start..range_x_end {
            let point_index = grid[index_fn(x, y)].checked_sub(1);
            if let Some(pi) = point_index {
                let sqr_dst = (candidate - points[pi]).length_squared();
                if sqr_dst <= r * r {
                    return false;
                }
            }
        }
    }

    if !is_in_polygon(candidate, polygon) {
        return false;
    }

    true
}
