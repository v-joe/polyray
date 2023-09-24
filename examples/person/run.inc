// Include file to make the Polyray person run
// Created by David Mason
// The transforms that are commented out can be added in to give
// the body a little smoother-looking run.

total_frames 25
start_frame 0
end_frame 24
outfile run

include "../colors.inc"
define personcolor shiny_red

// Steps is the number of cycles that will be rendered
// 2 actual steps is a "step"
define steps 1

define stepsize total_frames/steps
define time fmod(frame,stepsize)/stepsize

if (time<0.5) { define percent (time*2) }
          else { define percent (1-((time-0.5)*2)) }

define head_rotation transform { rotate <0,0,0> }
define neck_rotation transform { rotate <0,0,0> }
//define neck_rotation transform { rotate <-15,20-(40*percent),0> }
define torso_rotation transform { rotate <0,0,0> }
//define torso_rotation transform { rotate <0,-20+(40*percent),0> }
define abdomen_rotation transform { rotate <0,0,0> }
//define abdomen_rotation transform { rotate <0,0,-10+(20*percent)> }

define left_upper_rotation transform { rotate <-60+(105*percent),0,0> }
//define left_upper_rotation
//transform { rotate <-60+(105*percent),20-(40*percent),10-(20*percent)> }
define left_lower_rotation transform { rotate <-75+(60*percent),0,0> }
define left_hand_rotation transform { rotate <0,0,0> }
define left_finger_rotation transform { rotate <0,0,0> }

define right_upper_rotation transform { rotate <45+(-105*percent),0,0> }
//define right_upper_rotation
//transform { rotate <45+(-105*percent),20-(40*percent),10-(20*percent)> }
define right_lower_rotation transform { rotate <-15+(-60*percent),0,0> }
define right_hand_rotation transform { rotate <0,0,0> }
define right_finger_rotation transform { rotate <0,0,0> }

define left_thigh_rotation transform { rotate <15+(-95*percent),0,0> }
define left_calf_rotation transform { rotate <120+(-120*percent),0,0> }
define left_foot_rotation transform { rotate <0,0,0> }

define right_thigh_rotation transform { rotate <-80+(95*percent),0,0> }
define right_calf_rotation transform { rotate <120*percent,0,0> }
define right_foot_rotation transform { rotate <0,0,0> }