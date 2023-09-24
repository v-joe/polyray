define plane_blue
texture {
   surface {
      ambient steel_blue, 0.2
      diffuse steel_blue, 0.6
      specular white, 0.5
      reflection white, 0.3
      microfacet Reitz 20
      }
   }

define missile
object {
     object { cylinder <0, 0, 0>, <11, 0, 0>, 1 }
   + object { cone <11, 0, 0>, 1, <13, 0, 0>, 0 }
   + object { polygon 3, <0,  1,  0>, <0,  3,  0>, < 2,  1,  0> }
   + object { polygon 3, <0, -1,  0>, <0, -3,  0>, < 2, -1,  0> }
   + object { polygon 3, <0,  0,  1>, <0,  0,  3>, < 2,  0,  1> }
   + object { polygon 3, <0,  0, -1>, <0,  0, -3>, < 2,  0, -1> }
   + object { polygon 3, <9,  1,  0>, <9,  3,  0>, <11,  1,  0> }
   + object { polygon 3, <9, -1,  0>, <9, -3,  0>, <11, -1,  0> }
   + object { polygon 3, <9,  0,  1>, <9,  0,  3>, <11,  0,  1> }
   + object { polygon 3, <9,  0, -1>, <9,  0, -3>, <11,  0, -1> }
   shiny_red
   }

define left_wing
object {
     object { polygon 4, <10, 0, 0>, <30, 0,  0>, <18, 0, 20>, <10, 0, 20> }
   + object { polygon 3, <30, 0, 0>, <27, 0,  5>, <58, 0, 0> }
   }

define right_wing
object {
     object { polygon 4, <10, 0, 0>, <30, 0,  0>, <18, 0,-20>, <10, 0,-20> }
   + object { polygon 3, <30, 0, 0>, <27, 0, -5>, <58, 0, 0> }
   }

define left_flap
object { polygon 4, <1, 0, 0>, <6, 0, 0>, <3, 0, 8>, <-2, 0, 8> }

define right_flap
object { polygon 4, <1, 0, 0>, <6, 0, 0>, <3, 0,-8>, <-2, 0,-8> }

define engine
object {
     object { cylinder <0, 0, 0>, <30, 0, 0>, 2 }
   + object { cone <25, 0, 0>, 0, <30, 0, 0>, 2 }
   + object { cone <0, 0, 0>, 2, <-2, 0, 0>, 3 }
   }

define fuselage
object {
     object { cylinder <0, 0, 0>, <80, 0, 0>, 3 }
   + object { cone <80, 0, 0>, 3, <100, 0, 0>, 0 }
   scale <1, 1.5, 1>
   }

// Put together an aircraft
define fighter_plane
object {
     fuselage
   + engine { translate <0,-2, 3> }
   + engine { translate <0,-2,-3> }
   + left_wing  { translate <5, 1, 2> }
   + right_wing { translate <5, 1,-2> }
   + missile { translate <15, 1, 23> }
   + missile { translate <15, 1,-23> }
   + left_flap  { translate <0, 0, 3> }
   + right_flap { translate <0, 0,-3> }
   + left_flap  { rotate <-60, 0, 0> translate <0, 2, 2> }
   + right_flap { rotate < 60, 0, 0> translate <0, 2,-2> }
   plane_blue
   translate <-50, 0, 0>
   scale <0.06, 0.06, 0.06>
   }
