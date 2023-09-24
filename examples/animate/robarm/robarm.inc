// These first four values are normally all that needs to be changed
// in order to position the arm the way you want it.  (Of course
// calculation of the correct values can be tough.)
define pincher_angle         degrees(robot_arm_theta)
define wrist_rotation_angles <90, 0, degrees(robot_arm_beta)>

// Sizes of the various parts of the robot arm
define stand_height       3
define wrist1_radius      0.3
define wrist1_length      0.4
define wrist_radius       0.3
define wrist_length       0.4
define section1_length    4
define section2_length    4
define pincher_length     2
define pincher_width      pincher_length/7
define pincher_height     pincher_length/12
define finger_length      1.5
define finger_axis_length pincher_length/6.9
define finger_axis_radius pincher_length/6

// stand_orientation moves and orients the entire robot arm
define stand_orientation
transform {
   rotate <0, 0, 0> // Valid from <0, -180, 0> to <0, 180, 0>
   }

// arm_orientation moves the entire arm around the support
define arm_orientation
transform {
   rotate <0, 0, -degrees(robot_arm_beta)>
   rotate <0, degrees(robot_arm_alpha), 0>
   translate <0, stand_height, 0>
   }

define wrist1_orientation
transform {
   rotate wrist_rotation_angles
   translate <section2_length, 0, 0>
   }

// rotation of the entire wrist
define wrist_orientation
transform {
   rotate <0, 0, 0>    // Valid from <-45, 0, -180> to <45, 0, 180>
   }

// arm_extension defines how far the second arm piece extends
// from the first pieces
define arm_extension
transform {
   translate <section2_offset, 0, 0> // Valid from <1, 0, 0> to <4, 0, 0>
   }

define top_finger_orientation
transform {
   rotate <0, 0, pincher_angle>
   translate <pincher_length, 0, 0>
   }

define bottom_finger_orientation
transform {
   rotate <0, 0, -pincher_angle>
   translate <pincher_length, 0, 0>
   }

// The two pincher orientations determine if the hand is open or closed
define top_pincher_orientation
transform {
   rotate <0, 0, -pincher_angle>
   translate <0, 0.1, 0>
   }

define bottom_pincher_orientation
transform {
   rotate <0, 0, pincher_angle>
   translate <0, -0.1, 0>
   }

define wrist_base_axis
object {
   // Axis for the wrist
     object { cylinder <0, 0, -wrist_length>, <0, 0, wrist_length>, wrist_radius }
   + object { disc <0, 0, -wrist_length>, <0, 0, -1>, wrist_radius }
   + object { disc <0, 0, wrist_length>, <0, 0, 1>, wrist_radius }
   brass
   }

define finger
object {
     object {
      cylinder <0, 0, -finger_axis_length>,
               <0, 0,  finger_axis_length>,
               finger_axis_radius
      brass
      }
   + object {
      disc <0, 0, -finger_axis_length>, <0, 0, -1>, finger_axis_radius
      brass
      }
   + object {
      disc <0, 0, finger_axis_length>, <0, 0, 1>, finger_axis_radius
      brass
      }
   + object {
      box <0, -pincher_height, -pincher_width>, 
          <finger_length, pincher_height, pincher_width>
      matte_blue
      }
   }

define pincher
object {
   box <0, -pincher_height, -pincher_width>, 
       <pincher_length, pincher_height,  pincher_width>
   matte_blue
   }

define wrist_and_hand
object {
   wrist_base_axis
   + object {
          pincher
        + finger { top_finger_orientation }
        top_pincher_orientation
        }
   + object {
          pincher
        + finger { bottom_finger_orientation }
        bottom_pincher_orientation
        }
   wrist_orientation
   }

define complete_hand
object {
     object { cylinder <0, -wrist1_length, 0>, <0, wrist1_length, 0>,
                       wrist1_radius
              brass }
   + object { disc <0, -wrist1_length, 0>, <0, -1, 0>, wrist1_radius brass }
   + object { disc <0,  wrist1_length, 0>, <0,  1, 0>, wrist1_radius brass }
   + wrist_and_hand { translate <wrist1_radius, 0, 0> }
   wrist1_orientation
   }

define arm_section2
object {
   // Base part of the arm
     object { cylinder <0, 0, 0>, <section2_length, 0, 0>, 0.2 }
   + object { disc <section2_length, 0, 0>, <1, 0, 0>, 0.2 }
   + complete_hand
   arm_extension
   chrome
   }

define arm_section_sleeve
object {
     object { cylinder <-0.5, 0, 0>, <0.5, 0, 0>, 0.4 }
   + object { disc <0.5, 0, 0>, <-1, 0, 0>, 0.4 }
   + object { disc <0.5, 0, 0>, <1, 0, 0>, 0.4 }
   brass
   }

define arm_counterweight
object {
     object { cylinder <0, 0, 0>, <section1_length/4, 0, 0>, 0.5 }
   + object { disc <0, 0, 0>, <-1, 0, 0>, 0.5 }
   + object { disc <section1_length/4, 0, 0>, <1, 0, 0>, 0.5 }
   brass
   }

define arm_section1
object {
   // Base part of the arm
     object { cylinder <-section1_length/1.5, 0, 0>,
                       <section1_length/1.5, 0, 0>, 0.3 }
   + object { disc <section1_length/1.5, 0, 0>, <0, 0, 1>, 0.3 }
   + arm_section_sleeve { translate <section1_length/1.5, 0, 0> }
   + arm_counterweight { translate <-section1_length/1.3, 0, 0> }
   chrome
   }

define arm_base
object {
   // Axis for the arm
     object { cylinder <0, 0, -0.5>, <0, 0, 0.5>, 0.5 }
   + object { disc <0, 0, -0.5>, <0, 0, -1>, 0.5 }
   + object { disc <0, 0, 0.5>, <0, 0, 1>, 0.5 }
   brass
   }

define arm
object {
   arm_base + arm_section1 + arm_section2
   arm_orientation
   }

define stand
object {
     object { box <-2, 0, -2>, <2, 1, 2> brass }
   + object { cylinder <0, 1, 0>, <0, 2, 0>, 1 brass }
   + object { disc <0, 2, 0>, <0, 1, 0>, 1 brass }
   + object { cylinder <0, 2, 0>, <0, 3, 0>, 0.3 chrome }
   + object { disc <0, 3, 0>, <0, 1, 0>, 0.3 chrome }
   }

define robot_arm
object { stand { rotate <0, 30, 0> } + arm stand_orientation }
