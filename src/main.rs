mod sph;

use ggez::{event, GameError};
use ggez::graphics::{self, Color};
use ggez::{Context, GameResult};
// use ggez::glam::*;


struct MainState {
    sim: sph::Fluid,
}

impl MainState {
    fn new() -> GameResult<MainState> {
        let sim = sph::Fluid::new(500);
        Ok(MainState { sim })
    }
}

impl event::EventHandler<ggez::GameError> for MainState {
    fn update(&mut self, _ctx: &mut Context) -> GameResult {
        // Update the simulation
        self.sim.update();
        Ok(())
    }

    fn draw(&mut self, ctx: &mut Context) -> GameResult {
        let mut canvas = graphics::Canvas::from_frame(ctx, Color::from_rgb(255, 255, 255));

        // Draw the particles
        for i in self.sim.positions.iter() {
            
            let world_pos = graphics::Rect::new(i.x, i.y, 10.0, 10.0);
            
            let circle = graphics::Mesh::new_circle(
                ctx,
                graphics::DrawMode::fill(),
                *i,
                2.5,
                0.1,
                Color::from_rgb(0, 0, 255),
            )?;
            canvas.draw(&circle, graphics::DrawParam::default());
        }


        // Draw the boundary
        let boundary = graphics::Mesh::new_rectangle(
            ctx,
            graphics::DrawMode::stroke(1.0),
            graphics::Rect::new(0.0, 0.0, self.sim.bound_box.x, self.sim.bound_box.y),
            Color::from_rgb(255, 0, 0),
        )?;

        canvas.draw(&boundary, graphics::DrawParam::default());

        canvas.finish(ctx)?;
        Ok(())
    }
}

fn main() {
    let cb = ggez::ContextBuilder::new("fluid_sim", "gautam8404")
        .window_setup(ggez::conf::WindowSetup::default().title("Fluid Sim"))
        .window_mode(ggez::conf::WindowMode::default().dimensions(800.0, 600.0));
    
    let (ctx, event_loop) = cb.build().unwrap();
    let state = MainState::new().unwrap();
    event::run(ctx, event_loop, state)

}
