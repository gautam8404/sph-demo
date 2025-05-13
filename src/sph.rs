use glam::Vec2;
use std::f32::consts::PI;
use rayon::prelude::*;

const TIME_STEP: f32 = 0.0005;
// SMOOTHING RADIUS
const SR: f32 = 16.0;
const SR_2: f32 = SR * SR;
const LIQUID_DENSITY: f32 = 300.0; // kg/m^3
const GAS_CONSTANT: f32 = 2000.0;
const VISCOSITY: f32 = 200.0;
const GRAVITY: Vec2 = Vec2::new(0.0, 9.81); // m/s^2
const MASS: f32 = 2.5; // kg
const COLLISION_DAMP: f32 = 0.5; // damping factor for collision
const RADIUS: f32 = 5.0; // radius of the particles
const MIN_DISTANCE: f32 = 2.0 * RADIUS; // minimum distance between particles

#[derive(Debug)]
pub struct Fluid {
    pub positions: Vec<Vec2>,
    velocities: Vec<Vec2>,
    densities: Vec<f32>,
    pressures: Vec<f32>,
    forces: Vec<Vec2>,
    pub bound_box: Vec2,
    pub num_particles: usize,
}

// poly6 kernel function
// poly6 is used to calculate the density of a particle based on the distance to other particles
// as mentioned in the SPH paper
// r: distance between 2 particles
// h: smoothing radius
// 315/(64 * pi * h^9) * { (h^2 - r^2)^3 0 <= r <= h}, 0 otherwise
fn poly6(r: f32) -> f32 {
    let r2 = r * r;
    if r2 < SR_2 {
        return 315.0 / (64.0 * PI * SR.powi(9)) * (SR_2 - r2).powi(3);
    }

    0.0
}

// spiky kernel function
// spiky is used to calculate the pressure between 2 particles based on the distance to other particles
// spiky is used instead of poly 6 as poly6 approaches 0 as the distance approaches 0
// 15/(pi * h^6) * { (h - r)^3  0 <= r <= h}, 0 otherwise
// we use derivative of spiky kernel to get rate of change of pressure
// 15/(pi * h^6) * 3(h - r)^2 * - 1
// -45/(pi * h^6) * (h - r)^2
fn spiky(r: f32) -> f32 {
    if r < SR && r > 0.0 {
        return -45.0 / (PI * SR.powi(6)) * (SR - r).powi(2);
    }
    0.0
}

// viscosity kernel function
// viscosity is used to calculate the viscosity between 2 particles based on the distance to other particles
// 45/(pi * h^6) * { (h - r)  0 <= r <= h}, 0 otherwise
fn visc(r: f32) -> f32 {
    if r < SR && r > 0.0 {
        return 45.0 / (PI * SR.powi(6)) * (SR - r);
    }
    0.0
}

impl Fluid {
    pub fn new(num_particles: usize) -> Fluid {
        let mut positions = Vec::with_capacity(num_particles);

        for i in 0..num_particles {
            let x = (i % 20) as f32 * 12.0 + 100.0;
            let y = (i / 20) as f32 * 12.0 + 100.0;
            positions.push(Vec2::new(x, y));
        }

        let velocities = vec![Vec2::ZERO; num_particles];
        let densities = vec![0.0; num_particles];
        let pressures = vec![0.0; num_particles];
        let forces = vec![Vec2::ZERO; num_particles];
        let bound_box = Vec2::new(400.0, 400.0);
        let num_particles = num_particles;

        Fluid {
            positions,
            velocities,
            densities,
            pressures,
            forces,
            bound_box,
            num_particles,
        }
    }

    // for each particle calculate density based on distance to other particles which are within the smoothing radius
    // to calculate density we will use poly6 kernel function
    // pressure calculation is given in paper as p = k(ρ − ρ0),
    // where k is the stiffness constant, ρ is the density of the particle, and ρ0 is the rest density (liquid density)
    fn calculate_density_pressure(&mut self) {
        // for i in 0..self.num_particles {
        //     let mut density = 0.0;
        //     for j in 0..self.num_particles {
        //         let dist = (self.positions[i] - self.positions[j]).length();
        //         if dist < SR {
        //             density += MASS * poly6(dist);
        //         }
        //     }
        //     self.densities[i] = density;
        //     self.pressures[i] = GAS_CONSTANT * (density - LIQUID_DENSITY).max(0.0); // p = k(ρ − ρ0),
        // }
        self.densities
            .par_iter_mut()
            .zip(self.pressures.par_iter_mut())
            .enumerate()
            .for_each(|(i, (density, pressure))| {
                let mut d = 0.0;
                for j in 0..self.num_particles {
                    let dist = (self.positions[i] - self.positions[j]).length();
                    if dist < SR {
                        d += MASS * poly6(dist);
                    }
                }
                *density = d;
                *pressure = GAS_CONSTANT * (d - LIQUID_DENSITY).max(0.0); // p = k(ρ − ρ0),
            });
    }

    // for each particle calculate forces based on distance to other particles which are within the smoothing radius
    // fPressure is given in paper as follows
    //  -∑ M (pi + pj) / 2 * density_j * spiky(ri - rj, smoothing radius)
    // where M is mass, pi and pj are the pressures of the particles, and density_j is the density of the jth particle
    // fViscosity is given in paper as follows
    //  μ * ∑ M * (vi - vj) / 2 * density_j * visc(ri - rj, smoothing radius)
    // where μ is the viscosity, vi and vj are the velocities of the particles, and density_j is the density of the jth particle
    // gravity we can simply add to the force
    fn compute_forces(&mut self) {

        self.forces
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, force)|{
                let mut f_pressure = Vec2::ZERO;
                let mut f_viscosity = Vec2::ZERO;
        
                for j in 0..self.num_particles {
                    if i == j {
                        continue;
                    }
                    let r = self.positions[i] - self.positions[j];
                    let dist = r.length();
        
                    // (h - r)  0 <= r <= h
                    if dist < SR && dist > 0.0 {
                        let dir = r.normalize(); // direction of the force
                        // we use negative direction of the force as we want to push the particles away from each other
                        let f_p = MASS * (self.pressures[i] + self.pressures[j])
                            / (2.0 * self.densities[j])
                            * spiky(dist);
                        f_pressure += -dir * f_p;
        
                        let near_force = -100.0 * (SR - dist);
                        f_pressure += near_force * -dir;

                        if dist < MIN_DISTANCE {
                            let repulsion = 10.0 * -dir / dist;
                            f_pressure += repulsion;

                        }
        
                        f_viscosity += VISCOSITY * MASS * (self.velocities[j] - self.velocities[i])
                            / (self.densities[j])
                            * visc(dist);
        
                    }
                }
        
                // add gravity to the force
                let f_gravity = GRAVITY * MASS / self.densities[i];
                *force = f_pressure + f_viscosity + f_gravity;
            });

        // for i in 0..self.num_particles {
        //     let mut f_pressure = Vec2::ZERO;
        //     let mut f_viscosity = Vec2::ZERO;
        //     let repulsion = Vec2::ZERO;
        // 
        //     for j in 0..self.num_particles {
        //         if i == j {
        //             continue;
        //         }
        //         let r = self.positions[i] - self.positions[j];
        //         let dist = r.length();
        // 
        //         // (h - r)  0 <= r <= h
        //         if dist < SR && dist > 0.0 {
        //             let dir = r.normalize(); // direction of the force
        //             // we use negative direction of the force as we want to push the particles away from each other
        //             let f_p = MASS * (self.pressures[i] + self.pressures[j])
        //                 / (2.0 * self.densities[j])
        //                 * spiky(dist);
        //             f_pressure += -dir * f_p;
        // 
        //             let near_force = -100.0 * (SR - dist);
        //             f_pressure += near_force * -dir;
        // 
        //             if dist < MIN_DISTANCE {
        //                 let repulsion = 10.0 * -dir / dist;
        //                 f_pressure += repulsion;
        // 
        //             }
        // 
        //             f_viscosity += VISCOSITY * MASS * (self.velocities[j] - self.velocities[i])
        //                 / (self.densities[j])
        //                 * visc(dist);
        // 
        //             // let eps = 0.1;
        //             // let vel = eps * (self.velocities[j] - self.velocities[i]) * poly6(dist);
        //             // self.velocities[i] += vel;
        //             // if dist < MIN_DISTANCE {
        //             //     let overlap = MIN_DISTANCE - dist;
        //             //     let correction = 0.5 * overlap * dir;
        //             //     self.positions[i] += correction;
        //             //     self.positions[j] -= correction;
        //             // }
        //         }
        // 
        // 
        //     }
            // add gravity to the force
            // let f_gravity = GRAVITY * MASS / self.densities[i];
            // self.forces[i] = f_pressure + f_viscosity + f_gravity;
        // }
    }

    // update the velocity and position of the particles based on the forces
    fn integrate(&mut self) {
        for i in 0..self.num_particles {
            // update velocity based on the forces
            let density = self.densities[i].max(0.0001);
            self.velocities[i] += TIME_STEP * self.forces[i] / density;
            // update position based on the velocity
            self.positions[i] += self.velocities[i] * TIME_STEP;

            // check for boundary collisions
            if self.positions[i].x < 0.0 {
                self.positions[i].x = 0.0;
                self.velocities[i].x *= -COLLISION_DAMP
            }
            if self.positions[i].x > self.bound_box.x {
                self.positions[i].x = self.bound_box.x;
                self.velocities[i].x *= -COLLISION_DAMP;
            }
            if self.positions[i].y < 0.0 {
                self.positions[i].y = 0.0;
                self.velocities[i].y *= -COLLISION_DAMP;
            }
            if self.positions[i].y > self.bound_box.y {
                self.positions[i].y = self.bound_box.y;
                self.velocities[i].y *= -COLLISION_DAMP;
            }
        }
    }

    // update the fluid simulation
    pub fn update(&mut self) {
        println!("update");
        self.calculate_density_pressure();
        self.compute_forces();
        self.integrate();
    }
}
