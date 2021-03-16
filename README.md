
# physbox-rs
A physics sandbox written in Rust.
A toy project created while learning the language, graphics drawn using the `piston_window` crate.

![Orbits](https://github.com/For-The-Wolf/physbox-rs/blob/master/images/orbits.png)

**Figure 1**: Particle emittors and a particle orbiting a blackhole.

![Fountains](https://github.com/For-The-Wolf/physbox-rs/blob/master/images/fountains.png)

**Figure 2**: Emiitors with global gravity and containing walls enabled.

![Tidal Stripping](https://github.com/For-The-Wolf/physbox-rs/blob/master/images/tidal-stripping.png)

**Figure 3** A two black hole system, one stripping particles from the others orbit.


## Building
Simply `cd` into `physbox-rs` after cloning into the repo.

Then, build and run with `cargo run --release` (make sure to use the release flag to build with optimisation).

## Controls
* `Space` iterates through: Global gravity {Off, Down, Up}.
* `Enter` toggles solid walls around the window (on by default).
* `r` remvoes all entities from the simulation.
* `right-click` remove entities under the mouse pointer.
* `left-click` and drag, insert an entity depending on the current insert mode:
..* `p` Particle insert mode, creates a particle with the drawn initial velocity vector.
..* `e` Emmiter insert mode, creates a particle emitter emitting particles with velocities distributed about the drawn vector.
..* `b` BlackHole insert mode, creates a black hole with mass proportional to the drawn radius.
..* `o` Obstacle insert mode, inserts static collision detecting objects (not yet implelemted).

## Some examples
![Video 1](https://user-images.githubusercontent.com/74791897/111236386-a448ad80-85ea-11eb-8eb9-722e7bf8a8a5.mp4)

**Video 1**: Drawing two particle emitters and inserting a blackhole between them.

![Video 12](https://user-images.githubusercontent.com/74791897/111236521-f2f64780-85ea-11eb-8eea-c84d1b05a622.mp4)

**Video 2**: A two black hole system, one stripping mass from the other.

![Video 3](https://user-images.githubusercontent.com/74791897/111236554-086b7180-85eb-11eb-994d-23b3ebbe9851.mp4)

**Video 3**: Removing black holes from a system, then introducing global gravity.



