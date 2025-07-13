use std::{
    f32,
    ops::{Add, Mul, Sub},
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec2 {
    pub x: f32,
    pub y: f32,
}

impl Vec2 {
    pub fn new(x: f32, y: f32) -> Self {
        Self { x, y }
    }

    pub fn from_angle(angle: f32) -> Self {
        let (sin, cos) = f32::sin_cos(angle);
        Self { x: cos, y: sin }
    }

    pub fn dot(&self, other: &Self) -> f32 {
        self.x * other.x + self.y * other.y
    }

    pub fn length_squared(&self) -> f32 {
        self.dot(self)
    }

    pub fn normalize(&self) -> Self {
        let magnitude = self.length_squared().sqrt();
        if magnitude == 0. {
            Vec2 { x: 0., y: 0. }
        } else {
            Vec2 {
                x: self.x / magnitude,
                y: self.y / magnitude,
            }
        }
    }
}

impl Add<Vec2> for Vec2 {
    type Output = Vec2;
    #[inline]
    fn add(self, rhs: Vec2) -> Vec2 {
        Self {
            x: self.x.add(rhs.x),
            y: self.y.add(rhs.y),
        }
    }
}

impl Sub<Vec2> for Vec2 {
    type Output = Vec2;
    #[inline]
    fn sub(self, rhs: Vec2) -> Vec2 {
        Self {
            x: self.x.sub(rhs.x),
            y: self.y.sub(rhs.y),
        }
    }
}

impl Mul<Vec2> for Vec2 {
    type Output = Vec2;
    #[inline]
    fn mul(self, rhs: Vec2) -> Vec2 {
        Self {
            x: self.x.mul(rhs.x),
            y: self.y.mul(rhs.y),
        }
    }
}

impl Mul<f32> for Vec2 {
    type Output = Vec2;
    #[inline]
    fn mul(self, rhs: f32) -> Vec2 {
        Self {
            x: self.x.mul(rhs),
            y: self.y.mul(rhs),
        }
    }
}
