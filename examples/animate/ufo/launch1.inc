//
// This file builds the UFO and the launch pad
//

// Get the basic shape definition for a UFO
include "ufo.inc"

// Define the key positions of the UFO over the course of the animation
define start_pos1 [<0, -1, 0>, <0, -1, 0>, <0, 2, 0>, <0, 5, 0>]
define end_pos1   [<0, -1, 0>, <0, 2, 0>, <0, 5, 0>, <0, 5, 100>]
define pos0 start_pos1[segment]
define pos1 end_pos1[segment]

// Define the location and orientation of the UFO
define UFO_Location1
   start_pos1[segment] + increment * (end_pos1[segment] - start_pos1[segment])
define UFO_Orientation1
   transform {
      scale <1,1,1>/50        // Make it unit size
      rotate <0, 0, 0>      // Orient the UFO, then
      translate UFO_Location1 // move it into position
      }

// Change the orientations of the barrels for the ufo.
define Gun_Orientation1_1 transform { rotate <90, 30, 0> }

// Now make a UFO with a gun assembly hanging off the top and bottom
if (complexity == 1)
object {
   UFO_Body +
   object {
      // Mounting rod for the gun barrels.
        Rod { scale <1, Mounting_Rod_Length, 1> }
      // The twin gun barrels
      + Gun_Barrels { Gun_Orientation1_1 }
      translate Gun_Assembly_Location1
      Gun_Texture }
   UFO_Orientation1
   }
else
object {
   disc <0, 0, 0>, <0, 1, 0>, 50 
   UFO_Orientation1
   Gun_Texture
   }

// Launch pad for the UFO - this should be lots more complex...
define key_block object { box <4.5, 0, -0.2>, <5.5, 1.25, 0.2> steely_blue }
object {
     object { cylinder <0, -10, 0>, <0, 1, 0>, 1.5 }
   + object { cylinder <0, -1, 0>, <0, 1, 0>, 5 }
   + object { disc <0, 1, 0>, <0, 1, 0>, 1.5, 5 }
   + key_block
   + key_block { rotate <0, 60, 0> }
   + key_block { rotate <0,120, 0> }
   + key_block { rotate <0,180, 0> }
   + key_block { rotate <0,240, 0> }
   + key_block { rotate <0,300, 0> }
   matte_grey
   }
