use piston_window::{
    Context, DrawState, MouseCursorEvent, PressEvent, ReleaseEvent, RenderEvent, UpdateEvent,
};
//use std::path::Path;
use rand::random;
use std::collections::VecDeque;
use std::f64::consts::PI;
use std::time::SystemTime;

const WIDTH: f64 = 1280.0;
const HEIGHT: f64 = 720.0;
const GRAVITY: [f64; 2] = [0.0, 750.0];
const BIG_G: f64 = 670.0;
const PARTICLE_LIFETIME: f64 = 10.0;

#[derive(PartialEq)]
enum InsertState {
    Particle,
    Blackhole,
    Emitter,
    Static,
}

#[derive(PartialEq)]
enum ContainerState {
    Open,
    Closed,
}

#[derive(PartialEq)]
enum ButtonState {
    Pressed,
    Released,
}

#[derive(PartialEq)]
enum GlobalGravity {
    Up,
    Down,
    NoGravity,
}

struct Line {
    start: [f64; 2],
    end: [f64; 2],
}

struct Mouse {
    location: [f64; 2],
    prev_location: [f64; 2],
}
struct BlackHole {
    location: [f64; 2],
    radius: f64,
    mass: f64,
    colour: [f32; 4],
}
struct Circle {
    position: [f64; 2],
    velocity: [f64; 2],
    acceleration: [f64; 2],
    colour: [f32; 4],
    radius: f64,
    lifetime: f64,
    coeff_rest: f64,
    age: f64,
}

#[derive(Default)]
struct CircleParams {
    r: Option<f64>,
    p: Option<[f64; 2]>,
    v: Option<[f64; 2]>,
    a: Option<[f64; 2]>,
    c: Option<[f32; 4]>,
    l: Option<f64>,
    cr: Option<f64>,
}

struct Inflow {
    position: [f64; 2],
    flowrate: f64,
    velocity: [f64; 2],
    colour: [f32; 4],
    spawned: u64,
    shouldve_spawned: f64, //lifetime: f64
}

impl ContainerState {
    fn switch(&mut self) {
        match self {
            ContainerState::Open => *self = ContainerState::Closed,
            ContainerState::Closed => *self = ContainerState::Open,
        }
    }
}

impl BlackHole {
    fn new(location: [f64; 2], mass: f64, radius: f64) -> BlackHole {
        let colour = [0.1, 0.0, 0.2, 1.0];
        BlackHole {
            location,
            mass,
            radius,
            colour,
        }
    }
}

impl GlobalGravity {
    fn switch(&self) -> GlobalGravity {
        match self {
            GlobalGravity::NoGravity => GlobalGravity::Down,
            GlobalGravity::Down => GlobalGravity::Up,
            GlobalGravity::Up => GlobalGravity::NoGravity,
        }
    }
    fn num(&self) -> i32 {
        match self {
            GlobalGravity::NoGravity => 0,
            GlobalGravity::Down => 1,
            GlobalGravity::Up => -1,
        }
    }
}

impl Line {
    fn new() -> Line {
        Line {
            start: [0.0, 0.0],
            end: [0.0, 0.0],
        }
    }
    fn as_arr(&self) -> [f64; 2] {
        let [x1, y1] = self.start;
        let [x2, y2] = self.end;
        [x2 - x1, y2 - y1]
    }
    fn abs(&self) -> f64 {
        let [x1, y1] = self.start;
        let [x2, y2] = self.end;
        ((x2 - x1).powi(2) + (y2 - y1).powi(2)).sqrt()
    }
}

impl Inflow {
    fn new(position: [f64; 2], flowrate: f64, velocity: [f64; 2]) -> Inflow {
        Inflow {
            position,
            flowrate,
            velocity,
            spawned: 0,
            shouldve_spawned: 0.0,
            colour: random_colour(),
        }
    }
    fn spawn(&mut self, circles: &mut VecDeque<Circle>, dt: f64) {
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
    fn new(
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

    fn new_random(params: CircleParams) -> Circle {
        let radius = params.r.unwrap_or_else(|| random::<f64>() * 95.0 + 5.0);
        let position = params.p.unwrap_or_else(|| random_pos(radius));
        let velocity = params.v.unwrap_or_else(random_vel);
        let acceleration = params.a.unwrap_or_else(|| GRAVITY);
        let colour = params.c.unwrap_or_else(random_colour);
        let lifetime = params.l.unwrap_or_else(|| PARTICLE_LIFETIME);
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

    fn update_position(&mut self, dt: f64, container_state: &ContainerState, wall_thickness: &f64) {
        //let coeff_rest: f64 = 0.2 + random::<f64>() * 0.5;
        let coeff_rest: f64 = self.coeff_rest;
        self.age += dt;
        let [x, y] = self.position;
        let [vx, vy] = self.velocity;
        let [ax, ay] = self.acceleration;
        let mut new_x = x + vx * dt; // % WIDTH;
        let mut new_y = y + vy * dt; // % HEIGHT;
        let mut new_vx = vx + ax * dt;
        let mut new_vy = vy + ay * dt;
        if *container_state == ContainerState::Closed {
            if new_x - self.radius <= *wall_thickness {
                new_vx = -vx * coeff_rest;
                new_vy = vy * coeff_rest;
                new_x = self.radius + wall_thickness;
            } else if new_x + self.radius >= WIDTH - wall_thickness {
                new_vx = -vx * coeff_rest;
                new_vy = vy * coeff_rest;
                new_x = WIDTH - self.radius - wall_thickness;
            }
            if new_y - self.radius <= *wall_thickness {
                new_vy = -vy * coeff_rest;
                new_vx = vx * coeff_rest;
                new_y = self.radius + wall_thickness;
            } else if new_y + self.radius >= HEIGHT - wall_thickness {
                new_vy = -vy * coeff_rest;
                new_vx = vx * coeff_rest;
                new_y = HEIGHT - self.radius - wall_thickness;
            }
        }
        self.position = [new_x, new_y];
        self.velocity = [new_vx, new_vy];
    }
}

fn march_line(
    init_pos: [f64; 2],
    init_vect: [f64; 2],
    n_lines: u32,
    holes: &Vec<BlackHole>,
    g_direction: &GlobalGravity,
    container_state: &ContainerState,
    wall_thickness: &f64,
) -> Vec<Line> {
    let dt: f64 = 0.01;
    let mut tracker: Circle = Circle::new_random(CircleParams {
        p: Some(init_pos),
        v: Some(init_vect),
        r: Some(2.0),
        cr: Some(0.45),
        ..CircleParams::default()
    });
    let mut line_segments: Vec<Line> = Vec::new();
    for _ in 0..n_lines {
        let start: [f64; 2] = tracker.position;
        tracker.acceleration = [0.0, 0.0];
        for bh in holes {
            apply_gravity(bh, &mut tracker);
        }
        let [mut gx, mut gy] = GRAVITY;
        let [pax, pay] = tracker.acceleration;
        gy = gy * (g_direction.num() as f64);
        gx = gx * (g_direction.num() as f64);
        tracker.acceleration = [pax + gx, pay + gy];
        tracker.update_position(dt, container_state, wall_thickness);
        let end: [f64; 2] = tracker.position;
        line_segments.push(Line { start, end });
    }
    line_segments
}
fn _show_insert_mode(_insert_mode: InsertState) {}
fn draw_particles<G: piston_window::Graphics>(
    circs: &VecDeque<Circle>,
    d_state: &DrawState,
    c: Context,
    g: &mut G,
    container_state: &ContainerState,
) {
    for circ in circs {
        let [x, y] = circ.position;
        if container_state == &ContainerState::Closed
            || (x > -circ.radius)
                && (x < WIDTH + circ.radius)
                && (y > -circ.radius)
                && (y < HEIGHT + circ.radius)
        {
            let [c_r, c_g, c_b, _] = circ.colour;
            let ellipse = piston_window::Ellipse::new([
                c_r,
                c_g,
                c_b,
                1.0 - (circ.age / circ.lifetime) as f32,
            ])
            .resolution(7);
            ellipse.draw(
                [
                    x - circ.radius,
                    y - circ.radius,
                    circ.radius * 2.0,
                    circ.radius * 2.0,
                ],
                d_state,
                c.transform,
                g,
            );
        }
    }
}
fn draw_background<G: piston_window::Graphics>(
    container_state: &ContainerState,
    bg_colour: [f32; 4],
    wall_colour: [f32; 4],
    wall_thickness: f64,
    c: Context,
    g: &mut G,
) {
    if container_state == &ContainerState::Open {
        piston_window::clear(bg_colour, g);
    } else {
        piston_window::clear(wall_colour, g);
        piston_window::polygon(
            bg_colour,
            &[
                [WIDTH - wall_thickness, HEIGHT - wall_thickness],
                [wall_thickness, HEIGHT - wall_thickness],
                [wall_thickness, wall_thickness],
                [WIDTH - wall_thickness, wall_thickness],
            ],
            c.transform,
            g,
        );
    }
}
fn draw_indicator<G: piston_window::Graphics>(
    inflow_vector: &Line,
    insert_state: &InsertState,
    holes: &Vec<BlackHole>,
    g_direction: &GlobalGravity,
    container_state: &ContainerState,
    wall_thickness: &f64,
    bg_colour: [f32; 4],
    d_state: &DrawState,
    c: &Context,
    g: &mut G,
) {
    let [x1, y1] = inflow_vector.start;
    let [x2, y2] = inflow_vector.end;
    if insert_state == &InsertState::Particle || insert_state == &InsertState::Emitter {
        let segs = march_line(
            inflow_vector.start,
            inflow_vector.as_arr(),
            500,
            holes,
            g_direction,
            container_state,
            wall_thickness,
        );
        for line_seg in segs {
            let [seg_x1, seg_y1] = line_seg.start;
            let [seg_x2, seg_y2] = line_seg.end;
            piston_window::line(
                [0.0, 0.0, 0.0, 0.5],
                2.0,
                [seg_x1, seg_y1, seg_x2, seg_y2],
                c.transform,
                g,
            );
        }
        let inflow_indicator = piston_window::Line::new([0.7, 0.1, 0.0, 1.0], 2.0);
        inflow_indicator.draw_arrow([x1, y1, x2, y2], 5.0, d_state, c.transform, g);
        piston_window::ellipse(
            [0.7, 0.1, 0.0, 1.0],
            [x1 - 5.0, y1 - 5.0, 10.0, 10.0],
            c.transform,
            g,
        );
    } else if insert_state == &InsertState::Blackhole {
        let rad = inflow_vector.abs();
        piston_window::ellipse(
            [0.0, 0.0, 0.0, 0.4],
            [x1 - rad, y1 - rad, rad * 2.0, rad * 2.0],
            c.transform,
            g,
        );
        let inner_rad = f64::abs(rad - 2.0);
        piston_window::ellipse(
            bg_colour,
            [
                x1 - inner_rad,
                y1 - inner_rad,
                inner_rad * 2.0,
                inner_rad * 2.0,
            ],
            c.transform,
            g,
        );
    }
}

fn draw_black_holes<G: piston_window::Graphics>(
    holes: &mut Vec<BlackHole>,
    c: &Context,
    g: &mut G,
) {
    for black_hole in holes {
        let [bx, by] = black_hole.location;
        let br = black_hole.radius;
        piston_window::ellipse(
            black_hole.colour,
            [bx - br, by - br, br * 2.0, br * 2.0],
            c.transform,
            g,
        );
    }
}

fn draw_gravity_mode(g_direction: &GlobalGravity) -> ([[f64; 2]; 4]) {
    let x_len: f64 = 50.0;
    let y_len: f64 = 50.0;
    let x: f64 = 50.0;
    let y: f64 = 50.0;
    let points: [[f64; 2]; 4];
    match g_direction {
        GlobalGravity::Up => {
            points = [
                [x, y - (y_len / 2.0)],
                [x, y - (y_len / 2.0)],
                [x + (x_len / 2.0), y + (y_len / 2.0)],
                [x - (x_len / 2.0), y + (y_len / 2.0)],
            ];
        }
        GlobalGravity::NoGravity => {
            points = [
                [x - (x_len / 2.0), y + (y_len / 4.0)],
                [x + (x_len / 2.0), y + (y_len / 4.0)],
                [x + (x_len / 2.0), y - (y_len / 4.0)],
                [x - (x_len / 2.0), y - (y_len / 4.0)],
            ];
        }
        GlobalGravity::Down => {
            points = [
                [x - (x_len / 2.0), y - (y_len / 2.0)],
                [x + (x_len / 2.0), y - (y_len / 2.0)],
                [x, y + (y_len / 2.0)],
                [x, y + (y_len / 2.0)],
            ];
        }
    }
    points
}
fn update_particles(
    circles: &mut VecDeque<Circle>,
    holes: &Vec<BlackHole>,
    g_direction: &GlobalGravity,
    container_state: &ContainerState,
    wall_thickness: &f64,
    dt: f64,
) {
    for circ in circles {
        circ.acceleration = [0.0, 0.0];
        for bh in holes {
            apply_gravity(bh, circ);
        }
        let [mut gx, mut gy] = GRAVITY;
        let [pax, pay] = circ.acceleration;
        gy = gy * (g_direction.num() as f64);
        gx = gx * (g_direction.num() as f64);
        circ.acceleration = [pax + gx, pay + gy];
        circ.update_position(dt, container_state, wall_thickness)
    }
}
fn apply_gravity(bh: &BlackHole, particle: &mut Circle) {
    let [px, py] = particle.position;
    let [pax, pay] = particle.acceleration;
    let [bx, by] = bh.location;
    let distance = ((px - bx).powi(2) + (py - by).powi(2)).sqrt();
    let dx = (px - bx) / distance;
    let dy = (py - by) / distance;
    let force = -(BIG_G * bh.mass) / distance.powi(2);
    particle.acceleration = [pax + (force * dx), pay + (force * dy)];
}

fn random_pos(circ_rad: f64) -> [f64; 2] {
    let x = random::<f64>() * (WIDTH - (circ_rad * 2.0)) + circ_rad;
    let y = random::<f64>() * (HEIGHT - (circ_rad * 2.0)) + circ_rad;
    [x, y]
}

fn random_vel() -> [f64; 2] {
    let mag: f64 = 300.0;
    let ang = random::<f64>() * PI + PI;
    [mag * ang.cos(), mag * ang.sin()]
}
fn random_colour() -> [f32; 4] {
    [
        random::<f32>() * 0.7,
        random::<f32>() * 0.7,
        random::<f32>() * 0.7,
        1.0,
    ]
}

fn main() {
    let mut insert_state: InsertState = InsertState::Emitter;
    let wall_thickness: f64 = 2.0;
    let wall_colour: [f32; 4] = [0.9, 0.4, 0.2, 1.0];
    let mut inflow_vector: Line = Line::new();
    let mut container_state: ContainerState = ContainerState::Closed;
    let mut lmb: ButtonState = ButtonState::Released;
    let print_fps: bool = false;
    let mut g_direction: GlobalGravity = GlobalGravity::NoGravity;
    let mut mouse = Mouse {
        location: [0.0, 0.0],
        prev_location: [0.0, 0.0],
    };
    let bg_colour = [0.91, 0.9, 0.8, 1.0];
    let mut screen: piston_window::PistonWindow =
        piston_window::WindowSettings::new("physbox-rs", [WIDTH, HEIGHT])
            .exit_on_esc(true)
            .build()
            .unwrap();
    let mut events = screen.events;
    let mut circles: VecDeque<Circle> = VecDeque::new();
    let mut particles: VecDeque<Circle> = VecDeque::new();
    let mut spawners: Vec<Inflow> = Vec::new();
    let mut holes: Vec<BlackHole> = Vec::new();
    let mut then = SystemTime::now();
    let d_state = piston_window::DrawState::new_alpha();
    while let Some(e) = events.next(&mut screen) {
        if e.render_args().is_some() {
            let circs = &circles;
            screen.draw_2d(&e, |c, g, _| {
                draw_background(
                    &container_state,
                    bg_colour,
                    wall_colour,
                    wall_thickness,
                    c,
                    g,
                );
                if lmb == ButtonState::Pressed {
                    draw_indicator(
                        &inflow_vector,
                        &insert_state,
                        &holes,
                        &g_direction,
                        &container_state,
                        &wall_thickness,
                        bg_colour,
                        &d_state,
                        &c,
                        g,
                    );
                }
                draw_black_holes(&mut holes, &c, g);
                draw_particles(circs, &d_state, c, g, &container_state);
                draw_particles(&particles, &d_state, c, g, &container_state);
                for inflow in &mut spawners {
                    let [ix, iy] = inflow.position;
                    piston_window::polygon(
                        inflow.colour,
                        &[
                            [ix - 5.0, iy + 5.0],
                            [ix - 5.0, iy - 5.0],
                            [ix + 5.0, iy - 5.0],
                            [ix + 5.0, iy + 5.0],
                        ],
                        c.transform,
                        g,
                    );
                }
                let g_points = draw_gravity_mode(&g_direction);
                let g_colour = {
                    if g_direction.num() == 0 {
                        [1.0, 0.0, 0.0, 0.7]
                    } else {
                        [0.0, 1.0, 0.0, 0.7]
                    }
                };
                piston_window::polygon(g_colour, &g_points, c.transform, g);
            });
            if print_fps {
                let now = SystemTime::now();
                match now.duration_since(then) {
                    Ok(n) => println!(
                        "{} balls, {} fps.",
                        circles.len(),
                        (1e9 / (n.as_nanos() as f64))
                    ),
                    Err(_) => panic!("Oh no, time has failed us"),
                }
                then = now;
            }
        }
        if let Some(u) = e.update_args() {
            mouse.prev_location = mouse.location;
            for spawner in &mut spawners {
                spawner.spawn(&mut circles, u.dt);
            }
            loop {
                match circles.front() {
                    Some(last) if last.age > last.lifetime => circles.pop_front(),
                    _ => break,
                };
            }
            update_particles(
                &mut circles,
                &holes,
                &g_direction,
                &container_state,
                &wall_thickness,
                u.dt,
            );
            update_particles(
                &mut particles,
                &holes,
                &g_direction,
                &container_state,
                &wall_thickness,
                u.dt,
            );
        }
        if let Some(m) = e.mouse_cursor_args() {
            let cursor = &mut mouse;
            cursor.location = m;
            inflow_vector.end = m;
        }
        if let Some(p) = e.press_args() {
            match p {
                piston_window::Button::Keyboard(piston_window::Key::Space) => {
                    g_direction = g_direction.switch();
                }
                piston_window::Button::Keyboard(piston_window::Key::Return) => {
                    container_state.switch();
                }
                piston_window::Button::Keyboard(piston_window::Key::R) => {
                    circles = VecDeque::<Circle>::new();
                    spawners = Vec::new();
                    holes = Vec::new();
                }
                piston_window::Button::Keyboard(piston_window::Key::B) => {
                    insert_state = InsertState::Blackhole;
                }
                piston_window::Button::Keyboard(piston_window::Key::E) => {
                    insert_state = InsertState::Emitter;
                }
                piston_window::Button::Keyboard(piston_window::Key::S) => {
                    insert_state = InsertState::Static;
                }
                piston_window::Button::Keyboard(piston_window::Key::P) => {
                    insert_state = InsertState::Particle;
                }

                piston_window::Button::Mouse(piston_window::MouseButton::Left) => {
                    inflow_vector.start = mouse.location;
                    lmb = ButtonState::Pressed;
                    //let mag: f64 = 400.0;
                    //let ang = random::<f64>()*PI*2.0;
                    //let [vx,vy] = [mag*ang.cos(),mag*ang.sin()];
                    //spawners.push(Inflow::new([x, y], 500.0, [vx, vy]));
                }
                piston_window::Button::Mouse(piston_window::MouseButton::Right) => {
                    let [x, y] = mouse.location;
                    let mut to_remove: Vec<usize> = Vec::new();
                    for (n, spawner) in holes.iter().enumerate() {
                        let [s_x, s_y] = spawner.location;
                        let dist = ((x - s_x).powi(2) + (y - s_y).powi(2)).sqrt();
                        if dist < 5.0 {
                            to_remove.push(n);
                        }
                    }
                    for n in to_remove.iter().rev() {
                        holes.remove(*n);
                    }
                }
                _ => (),
            }
        }
        if let Some(r) = e.release_args() {
            match r {
                piston_window::Button::Mouse(piston_window::MouseButton::Left) => {
                    inflow_vector.end = mouse.location;
                    lmb = ButtonState::Released;
                    match insert_state {
                        InsertState::Emitter => {
                            spawners.push(Inflow::new(
                                inflow_vector.start,
                                250.0,
                                inflow_vector.as_arr(),
                            ));
                        }
                        InsertState::Blackhole => {
                            holes.push(BlackHole::new(
                                inflow_vector.start,
                                inflow_vector.abs() * 100.0,
                                5.0,
                            ));
                        }
                        InsertState::Particle => {
                            particles.push_back(Circle::new(
                                inflow_vector.start,
                                inflow_vector.as_arr(),
                                [0.0, 0.0],
                                10.0,
                                random_colour(),
                                0.45,
                                100_000.0,
                            ));
                        }
                        _ => (),
                    }
                }
                _ => (),
            }
        }
    }
}
