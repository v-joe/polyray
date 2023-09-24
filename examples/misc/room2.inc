// This is a scene used in room1.pi.

// Make the test scene - cylinders arranged in a hexagon on a checkered plane.
// Inside the cylinders will be a reflective sphere.
define bot_vec <6, -2, 0>
define top_vec <6, 8, 0>
define cyl_obj
object {
     object { cylinder bot_vec, top_vec, 1 }
   + object { disc top_vec, top_vec - bot_vec, 1 uv_steps 8, 4 }
   }
cyl_obj { rotate <0,   0, 0> shiny_red }
cyl_obj { rotate <0,  60, 0> shiny_yellow }
cyl_obj { rotate <0, 120, 0> shiny_green }
cyl_obj { rotate <0, 180, 0> shiny_cyan }
cyl_obj { rotate <0, 240, 0> shiny_blue }
cyl_obj { rotate <0, 300, 0> shiny_magenta }

object {
   disc <0, -2.001, 0>, <0, 1, 0>, 50
   texture { checker matte_white, matte_black }
   }
