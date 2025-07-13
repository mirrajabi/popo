use std::cmp::Ordering;

use crate::vectors::Vec2;

pub type Polygon = Vec<Vec2>;

/// Offsets a polygon by a given distance `value`.
/// This function mutates the input polygon to ensure it is in counter-clockwise order before offsetting.
/// A positive `value` offsets it outwards [expand], while a negative `value` offsets the polygon inwards [shrink].
pub fn offset(mut polygon: Polygon, value: f32) -> Polygon {
    polygon = sort_ccw(polygon);

    let n_verts = polygon.len();
    let mut offset_polygon = Vec::with_capacity(n_verts);
    for i in 0..n_verts {
        let prev = polygon[(i + n_verts - 1) % n_verts];
        let pi = polygon[i];
        let next = polygon[(i + 1) % n_verts];

        let l1 = (pi - prev).normalize();
        let l2 = (next - pi).normalize();

        let n1 = Vec2 { x: l1.y, y: -l1.x };
        let n2 = Vec2 { x: l2.y, y: -l2.x };

        let bisector = (n1 + n2).normalize();
        // Angle between edges
        let dot = (l1 * -1.).dot(&l2).max(-1.).min(1.);
        let theta = f32::acos(dot);
        let offset_length = value / f32::sin(theta / 2.);
        let new_pi = pi + bisector * offset_length;

        offset_polygon.push(new_pi);
    }

    offset_polygon
}

/// Sorts a given polygon in CCW order.  
/// If the input is already in CCW order or is colinear, it will be returned as is. Otherwise it will be cloned and reversed.
fn sort_ccw(polygon: Polygon) -> Polygon {
    if polygon.len() < 3 {
        return polygon;
    }

    let winding = find_winding_order(&polygon);

    // In case of colinear, there's no space to fit any circles anyway, so why bother! Just return the polygon
    if winding.is_ge() {
        polygon
    } else {
        let mut out = polygon.clone();
        out.reverse();
        out
    }
}

/// Finds the winding order (CW/CCW) of the polygon by calculating the signed area of a polygon and comparing it to zero.
/// If the output:
/// - `is_gt` 0, the polygon is in CCW order,
/// - `is_eq` 0, the polygon is Colinear (ordering can't be determined; Inwards and outwards have no meaning),
/// - `is_lt` 0, the polygon is in CW order
fn find_winding_order(polygon: &Polygon) -> Ordering {
    find_signed_area(polygon).total_cmp(&0.)
}

/// Finds the signed area of a polygon using the [Shoelace formula (Trapezoid formula)](https://en.wikipedia.org/wiki/Shoelace_formula)
fn find_signed_area(polygon: &Polygon) -> f32 {
    let n_verts = polygon.len();

    let mut area = 0.;
    for i in 0..n_verts {
        let current = polygon[i];
        let next = polygon[(i + 1) % n_verts];

        area += (current.x * next.y) - (next.x * current.y);
    }

    area / 2.
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vectors::Vec2;

    #[test]
    fn test_sort_ccw_two_points_polygon() {
        let line = vec![Vec2 { x: 0.0, y: 0.0 }, Vec2 { x: 1.0, y: 0.0 }];
        let sorted = sort_ccw(line.clone());
        assert_eq!(sorted, line);
    }

    #[test]
    fn test_sort_ccw_ccw_polygon() {
        // Triangle below has vertices in counterclockwise order.
        // Points (0, 0), (1, 0), (0, 1) yield a positive area.
        let ccw_polygon = vec![
            Vec2 { x: 0.0, y: 0.0 },
            Vec2 { x: 1.0, y: 0.0 },
            Vec2 { x: 0.0, y: 1.0 },
        ];
        let sorted = sort_ccw(ccw_polygon.clone());
        // Since the polygon is already counter-clockwise, its order remains unchanged.
        assert_eq!(sorted, ccw_polygon);
    }

    #[test]
    fn test_sort_ccw_cw_polygon() {
        // Triangle below has vertices in clockwise order.
        // Points (0, 1), (1, 0), (0, 0) yield a negative area.
        let cw_polygon = vec![
            Vec2 { x: 0.0, y: 1.0 },
            Vec2 { x: 1.0, y: 0.0 },
            Vec2 { x: 0.0, y: 0.0 },
        ];
        let sorted = sort_ccw(cw_polygon.clone());
        // When winding is non-negative (CCW or colinear), we reverse the vertices.
        let mut expected = cw_polygon;
        expected.reverse();
        assert_eq!(sorted, expected);
    }

    #[test]
    fn test_sort_ccw_colinear_polygon() {
        // Polygon with an area equal to 0
        let colinear = vec![
            Vec2 { x: 0.0, y: 0.0 },
            Vec2 { x: 1.0, y: 1.0 },
            Vec2 { x: 2.0, y: 2.0 },
        ];
        let sorted = sort_ccw(colinear.clone());
        // Colinear points are treated as not needing a reversal.
        assert_eq!(sorted, colinear);
    }

    #[test]
    fn test_find_signed_area() {
        let polygon = vec![
            Vec2 { x: 0.0, y: 0.0 },
            Vec2 { x: 1.0, y: 0.0 },
            Vec2 { x: 0.0, y: 1.0 },
        ];
        let area = find_signed_area(&polygon);
        // The area of the triangle is 0.5
        assert_eq!(area, 0.5);
    }

    #[test]
    fn test_find_signed_area_negative() {
        let polygon = vec![
            Vec2 { x: 0.0, y: 1.0 },
            Vec2 { x: 1.0, y: 0.0 },
            Vec2 { x: 0.0, y: 0.0 },
        ];
        let area = find_signed_area(&polygon);
        // The area of the triangle is -0.5
        assert_eq!(area, -0.5);
    }

    #[test]
    fn test_find_signed_area_colinear() {
        let polygon = vec![
            Vec2 { x: 0.0, y: 0.0 },
            Vec2 { x: 1.0, y: 1.0 },
            Vec2 { x: 2.0, y: 2.0 },
        ];
        let area = find_signed_area(&polygon);
        // The area of colinear points is 0
        assert_eq!(area, 0.0);
    }

    #[test]
    fn test_offset_polygon() {
        let polygon = vec![
            Vec2 { x: 0.0, y: 0.0 },
            Vec2 { x: 1.0, y: 0.0 },
            Vec2 { x: 0.0, y: 1.0 },
        ];
        let offset_value = 0.1;
        let offset_polygon = offset(polygon, offset_value);
        // The offset polygon should have vertices shifted outward by the offset value.
        assert_eq!(offset_polygon.len(), 3);
        assert!(offset_polygon[0].x < 0. && offset_polygon[0].y < 0.);
        assert!(offset_polygon[1].x > 1. && offset_polygon[1].y < 0.);
        assert!(offset_polygon[2].x < 0. && offset_polygon[2].y > 1.);
    }

    #[test]
    fn test_offset_polygon_negative() {
        let polygon = vec![
            Vec2 { x: 0.0, y: 0.0 },
            Vec2 { x: 1.0, y: 0.0 },
            Vec2 { x: 0.0, y: 1.0 },
        ];
        let offset_value = -0.1;
        let offset_polygon = offset(polygon, offset_value);
        // The offset polygon should have vertices shifted inward by the offset value.
        assert_eq!(offset_polygon.len(), 3);
        assert!(offset_polygon[0].x > 0. && offset_polygon[0].y > 0.);
        assert!(offset_polygon[1].x < 1. && offset_polygon[1].y > 0.);
        assert!(offset_polygon[2].x > 0. && offset_polygon[2].y < 1.);
    }
}
