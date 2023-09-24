viewpoint {
   from <0, 10, -14>
   at <0,0,0>
   up <0,1,0>
   angle 45
   resolution 160, 120
   aspect 4/3
   }

background <0, 0, 0>
light <-50, 10,  50>
light < 10, 20, 100>

include "../colors.inc"

define blue_ripple
texture {
   special shiny {
      color <0.4, 0.4, 1.0>
      normal N + ripple_fn
      reflection 0.5
      }
   }

object {
   disc <0, 0, 0>, <0, 1, 0>, 300
   blue_ripple { scale <3, 5, 3> }
   }

