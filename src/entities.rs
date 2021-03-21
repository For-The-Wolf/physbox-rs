use super::consts;
use rand::random;
use std::collections::VecDeque;
use std::f64::consts::PI;

pub struct Line {
    pub start: [f64; 2],
    pub end: [f64; 2],
}
pub struct BlackHole {
    pub location: [f64; 2],
    pub radius: f64,
    pub mass: f64,
    pub colour: [f32; 4],
}
pub struct Circle {
    pub position: [f64; 2],
    pub velocity: [f64; 2],
    pub acceleration: [f64; 2],
    pub colour: [f32; 4],
    pub radius: f64,
    pub lifetime: f64,
    pub coeff_rest: f64,
    pub age: f64,
}

#[derive(Default)]
pub struct CircleParams {
    pub r: Option<f64>,
    pub p: Option<[f64; 2]>,
    pub v: Option<[f64; 2]>,
    pub a: Option<[f64; 2]>,
    pub c: Option<[f32; 4]>,
    pub l: Option<f64>,
    pub cr: Option<f64>,
}

pub struct Inflow {
    pub position: [f64; 2],
    pub flowrate: f64,
    pub velocity: [f64; 2],
    pub colour: [f32; 4],
    spawned: u64,
    shouldve_spawned: f64, //lifetime: f64
}
#[derive(PartialEq)]
pub enum ContainerState {
    Open,
    Closed,
}
impl BlackHole {
    pub fn new(location: [f64; 2], mass: f64, radius: f64) -> BlackHole {
        let colour = [0.1, 0.0, 0.2, 1.0];
        BlackHole {
            location,
            mass,
            radius,
            colour,
        }
    }
}
impl ContainerState {
    pub fn switch(&mut self) {
        match self {
            ContainerState::Open => *self = ContainerState::Closed,
            ContainerState::Closed => *self = ContainerState::Open,
        }
    }
}

impl Line {
    pub fn new() -> Line {
        Line {
            start: [0.0, 0.0],
            end: [0.0, 0.0],
        }
    }
    pub fn as_arr(&self) -> [f64; 2] {
        let [x1, y1] = self.start;
        let [x2, y2] = self.end;
        [x2 - x1, y2 - y1]
    }
    pub fn abs(&self) -> f64 {
        let [x1, y1] = self.start;
        let [x2, y2] = self.end;
        ((x2 - x1).powi(2) + (y2 - y1).powi(2)).sqrt()
    }
}

impl Inflow {
    pub fn new(position: [f64; 2], flowrate: f64, velocity: [f64; 2]) -> Inflow {
        Inflow {
            position,
            flowrate,
            velocity,
            spawned: 0,
            shouldve_spawned: 0.0,
            colour: random_colour(),
        }
    }
    pub fn spawn(&mut self, circles: &mut VecDeque<Circle>, dt: f64) {
        self.shouldve_spawned += self.flowrate * dt;
        let [vx, vy] = self.velocity;
        let mag = (vx.powi(2) + vy.powi(2)).sqrt();
        let dir = vy.atan2(vx);
        let [r, g, b, _] = self.colour;
        while self.spawned as f64 <= self.shouldve_spawned {
            let c_mag = mag + (random::<f64>() - 0.5) * mag * 0.05;
            let c_dir = dir + (random::<f64>() - 0.5) * PI / 6.0;
            let c_vy = c_mag * c_dir.sin();
            let c_vx = c_mag * c_dir.cos();
            circles.push_back(Circle::new_random(CircleParams {
                r: Some(6.0),
                p: Some(self.position),
                v: Some([c_vx, c_vy]),
                c: Some([r, g, b, 0.7]),
                ..CircleParams::default()
            }));
            self.spawned += 1;
        }
    }
}
impl Circle {
    pub fn new(
        position: [f64; 2],
        velocity: [f64; 2],
        acceleration: [f64; 2],
        radius: f64,
        colour: [f32; 4],
        coeff_rest: f64,
        lifetime: f64,
    ) -> Circle {
        Circle {
            position,
            velocity,
            acceleration,
            radius,
            colour,
            coeff_rest,
            lifetime,
            age: 0.0,
        }
    }

    pub fn new_random(params: CircleParams) -> Circle {
        let radius = params.r.unwrap_or_else(|| random::<f64>() * 95.0 + 5.0);
        let position = params.p.unwrap_or_else(|| random_pos(radius));
        let velocity = params.v.unwrap_or_else(random_vel);
        let acceleration = params.a.unwrap_or_else(|| consts::GRAVITY);
        let colour = params.c.unwrap_or_else(random_colour);
        let lifetime = params.l.unwrap_or_else(|| consts::PARTICLE_LIFETIME);
        let coeff_rest = params.cr.unwrap_or_else(|| random::<f64>() * 0.5 + 0.2);
        Circle::new(
            position,
            velocity,
            acceleration,
            radius,
            colour,
            coeff_rest,
            lifetime,
        )
    }

    pub fn update_position(
        &mut self,
        dt: f64,
        container_state: &ContainerState,
        wall_thickness: &f64,
    ) {
        //let coeff_rest: f64 = 0.2 + random::<f64>() * 0.5;
        let coeff_rest: f64 = self.coeff_rest;
        self.age += dt;
        let [x, y] = self.position;
        let [vx, vy] = self.velocity;
        let [ax, ay] = self.acceleration;
        let mut new_x = x + vx * dt; // % consts::WIDTH;
        let mut new_y = y + vy * dt; // % consts::HEIGHT;
        let mut new_vx = vx + ax * dt;
        let mut new_vy = vy + ay * dt;
        if *container_state == ContainerState::Closed {
            if new_x - self.radius <= *wall_thickness {
                new_vx = -vx * coeff_rest;
                new_vy = vy * coeff_rest;
                new_x = self.radius + wall_thickness;
            } else if new_x + self.radius >= consts::WIDTH - wall_thickness {
                new_vx = -vx * coeff_rest;
                new_vy = vy * coeff_rest;
                new_x = consts::WIDTH - self.radius - wall_thickness;
            }
            if new_y - self.radius <= *wall_thickness {
                new_vy = -vy * coeff_rest;
                new_vx = vx * coeff_rest;
                new_y = self.radius + wall_thickness;
            } else if new_y + self.radius >= consts::HEIGHT - wall_thickness {
                new_vy = -vy * coeff_rest;
                new_vx = vx * coeff_rest;
                new_y = consts::HEIGHT - self.radius - wall_thickness;
            }
        }
        self.position = [new_x, new_y];
        self.velocity = [new_vx, new_vy];
    }
}
pub fn apply_gravity(bh: &BlackHole, particle: &mut Circle) {
    let [px, py] = particle.position;
    let [pax, pay] = particle.acceleration;
    let [bx, by] = bh.location;
    let distance = ((px - bx).powi(2) + (py - by).powi(2)).sqrt();
    let dx = (px - bx) / distance;
    let dy = (py - by) / distance;
    let force = -(consts::BIG_G * bh.mass) / distance.powi(2);
    particle.acceleration = [pax + (force * dx), pay + (force * dy)];
}

fn random_pos(circ_rad: f64) -> [f64; 2] {
    let x = random::<f64>() * (consts::WIDTH - (circ_rad * 2.0)) + circ_rad;
    let y = random::<f64>() * (consts::HEIGHT - (circ_rad * 2.0)) + circ_rad;
    [x, y]
}

fn random_vel() -> [f64; 2] {
    let mag: f64 = 300.0;
    let ang = random::<f64>() * PI + PI;
    [mag * ang.cos(), mag * ang.sin()]
}
pub fn random_colour() -> [f32; 4] {
    [
        random::<f32>() * 0.7,
        random::<f32>() * 0.7,
        random::<f32>() * 0.7,
        1.0,
    ]
}
