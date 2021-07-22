use piston_window::{
    Context, DrawState, MouseCursorEvent, PressEvent, ReleaseEvent, RenderEvent, UpdateEvent,
};
//use std::path::Path;
use std::collections::VecDeque;
use std::time::SystemTime;
mod consts;
mod entities;

#[derive(PartialEq)]
enum InsertState {
    Particle,
    Blackhole,
    Emitter,
    Static,
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

struct Mouse {
    location: [f64; 2],
    prev_location: [f64; 2],
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

fn march_line(
    init_pos: [f64; 2],
    init_vect: [f64; 2],
    n_lines: u32,
    holes: &Vec<entities::BlackHole>,
    g_direction: &GlobalGravity,
    container_state: &entities::ContainerState,
    wall_thickness: &f64,
) -> Vec<entities::Line> {
    let dt: f64 = 0.01;
    let mut tracker: entities::Circle = entities::Circle::new_random(entities::CircleParams {
        p: Some(init_pos),
        v: Some(init_vect),
        r: Some(2.0),
        cr: Some(0.45),
        ..entities::CircleParams::default()
    });
    let mut line_segments: Vec<entities::Line> = Vec::new();
    for _ in 0..n_lines {
        let start: [f64; 2] = tracker.position;
        tracker.acceleration = [0.0, 0.0];
        for bh in holes {
            entities::apply_gravity(bh, &mut tracker);
        }
        let [mut gx, mut gy] = consts::GRAVITY;
        let [pax, pay] = tracker.acceleration;
        gy = gy * (g_direction.num() as f64);
        gx = gx * (g_direction.num() as f64);
        tracker.acceleration = [pax + gx, pay + gy];
        tracker.update_position(dt, container_state, wall_thickness);
        let end: [f64; 2] = tracker.position;
        line_segments.push(entities::Line { start, end });
    }
    line_segments
}

fn _show_insert_mode(_insert_mode: InsertState) {
    todo!();
}

fn draw_particles<G: piston_window::Graphics>(
    circs: &VecDeque<entities::Circle>,
    d_state: &DrawState,
    c: Context,
    g: &mut G,
    container_state: &entities::ContainerState,
) {
    for circ in circs {
        let [x, y] = circ.position;
        if container_state == &entities::ContainerState::Closed
            || (x > -circ.radius)
                && (x < consts::WIDTH + circ.radius)
                && (y > -circ.radius)
                && (y < consts::HEIGHT + circ.radius)
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
    container_state: &entities::ContainerState,
    bg_colour: [f32; 4],
    wall_colour: [f32; 4],
    wall_thickness: f64,
    c: Context,
    g: &mut G,
) {
    if container_state == &entities::ContainerState::Open {
        piston_window::clear(bg_colour, g);
    } else {
        piston_window::clear(wall_colour, g);
        piston_window::polygon(
            bg_colour,
            &[
                [
                    consts::WIDTH - wall_thickness,
                    consts::HEIGHT - wall_thickness,
                ],
                [wall_thickness, consts::HEIGHT - wall_thickness],
                [wall_thickness, wall_thickness],
                [consts::WIDTH - wall_thickness, wall_thickness],
            ],
            c.transform,
            g,
        );
    }
}

fn draw_indicator<G: piston_window::Graphics>(
    inflow_vector: &entities::Line,
    insert_state: &InsertState,
    holes: &Vec<entities::BlackHole>,
    g_direction: &GlobalGravity,
    container_state: &entities::ContainerState,
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
    holes: &mut Vec<entities::BlackHole>,
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
    circles: &mut VecDeque<entities::Circle>,
    holes: &Vec<entities::BlackHole>,
    g_direction: &GlobalGravity,
    container_state: &entities::ContainerState,
    wall_thickness: &f64,
    dt: f64,
) {
    for circ in circles {
        circ.acceleration = [0.0, 0.0];
        for bh in holes {
            entities::apply_gravity(bh, circ);
        }
        let [mut gx, mut gy] = consts::GRAVITY;
        let [pax, pay] = circ.acceleration;
        gy = gy * (g_direction.num() as f64);
        gx = gx * (g_direction.num() as f64);
        circ.acceleration = [pax + gx, pay + gy];
        circ.update_position(dt, container_state, wall_thickness)
    }
}

fn main() {
    let mut insert_state: InsertState = InsertState::Emitter;
    let wall_thickness: f64 = 2.0;
    let wall_colour: [f32; 4] = [0.9, 0.4, 0.2, 1.0];
    let mut inflow_vector: entities::Line = entities::Line::new();
    let mut container_state: entities::ContainerState = entities::ContainerState::Closed;
    let mut lmb: ButtonState = ButtonState::Released;
    let print_fps: bool = false;
    let mut g_direction: GlobalGravity = GlobalGravity::NoGravity;
    let mut mouse = Mouse {
        location: [0.0, 0.0],
        prev_location: [0.0, 0.0],
    };
    let bg_colour = [0.91, 0.9, 0.8, 1.0];
    let mut screen: piston_window::PistonWindow =
        piston_window::WindowSettings::new("physbox-rs", [consts::WIDTH, consts::HEIGHT])
            .exit_on_esc(true)
            .build()
            .unwrap();
    let mut events = screen.events;
    let mut circles: VecDeque<entities::Circle> = VecDeque::new();
    let mut particles: VecDeque<entities::Circle> = VecDeque::new();
    let mut spawners: Vec<entities::Inflow> = Vec::new();
    let mut holes: Vec<entities::BlackHole> = Vec::new();
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
                    circles = VecDeque::<entities::Circle>::new();
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
                    //spawners.push(entities::Inflow::new([x, y], 500.0, [vx, vy]));
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
                            spawners.push(entities::Inflow::new(
                                inflow_vector.start,
                                250.0,
                                inflow_vector.as_arr(),
                            ));
                        }
                        InsertState::Blackhole => {
                            holes.push(entities::BlackHole::new(
                                inflow_vector.start,
                                inflow_vector.abs() * 100.0,
                                5.0,
                            ));
                        }
                        InsertState::Particle => {
                            particles.push_back(entities::Circle::new(
                                inflow_vector.start,
                                inflow_vector.as_arr(),
                                [0.0, 0.0],
                                10.0,
                                entities::random_colour(),
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
