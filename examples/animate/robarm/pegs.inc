define peg_center_texture light_wood { rotate <75, 0, 0> }
define peg_base_texture light_wood { rotate <30, 30, 0> }

define peg_length      2.5
define peg_width       0.4
define peg_base_width  1.5
define peg_base_height 0.4
define peg_offset      4
define peg_distance    8
define peg_center
object {
    object { sphere <0, peg_length, 0>, peg_width }
  + object { cylinder <0, 0, 0>, <0, peg_length, 0>, peg_width }
  }
define peg_base
object {
    object { cylinder <0, 0, 0>, <0, peg_base_height, 0>, peg_base_width }
  + object { disc <0, peg_base_height, 0>, <0, 1, 0>, peg_base_width }
  }
define pegs
object {
   // Peg #1
   object {
        peg_center {
          peg_center_texture
          }
      + peg_base {
          peg_base_texture
          }
      translate <0, 0, -peg_offset>
      }
   // Peg #2
 + object {
        peg_center {
          peg_center_texture {
            translate <0, 0, 20>
            rotate <10, 0, -10>
            }
          }
      + peg_base {
          peg_base_texture {
             translate <0, 0, 20>
             rotate <5, -50, 0>
             }
          }
     }
   // Peg #3
 + object {
        peg_center {
          peg_center_texture {
            translate <0, 0, 40>
            rotate <0, 30, -10>
            }
          }
      + peg_base {
          peg_base_texture {
             translate <0, 0, 40>
             rotate <15, 0, -20>
             }
          }
     translate <0, 0, peg_offset>
     }
  }
