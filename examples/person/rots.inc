// Polyray include file for person.pi
// This file contains the rotations for all the joints of the body,
// and color of the figure.
// All rotations follow the left hand rule.
// +x is to the right, +z is to the front, and +y is up.
// Created by Will Wagner.
// Thanks to Alexander Enzmann for his help.

define head_rotation		// <-45/45-80/80,-10/10,>
transform {
	rotate <0,0,0>
	}
define neck_rotation		// <-30/30,-20/20,-10/10>
transform {
	rotate <0,0,0>
	}
define torso_rotation		// <-10/100,-15/15,-25/25>
transform {
	rotate <0,0,0>
	}
define abdomen_rotation		// <0,0,-5/5>
transform {
	rotate <0,0,0>
	}

// Left arm rotations
define left_upper_rotation		// <-180/45,90/-45,45/-180>
transform {
	rotate <0,0,0>
	}
define left_lower_rotation		// <-165/0,0,0>
transform {
	rotate <-10,0,0>
	}
define left_hand_rotation		// <-20/20,90/-180,90/-20>
transform {
	rotate <0,30,0>
	}
define left_finger_rotation		// <0,0,90/0>
transform {
	rotate <0,0,45>
	}

// Right arm rotations
define right_upper_rotation		// <-180/45,45/-90,180/-45>
transform {
	rotate <0,0,0>
	}
define right_lower_rotation		// <-165/0,0,0>
transform {
	rotate <-10,0,0>
	}
define right_hand_rotation		// <-20/20,180/-90,20/-90>
transform {
	rotate <0,-30,0>
	}
define right_finger_rotation		// <0,0,0/-90>
transform {
	rotate <0,0,-45>
	}

// Left leg rotations
define left_thigh_rotation		// <-90/15,45/0,0/-45>
transform {
	rotate <0,0,0>
	}
define left_calf_rotation		// <0/160,0,0>
transform {
	rotate <0,0,0>
	}
define left_foot_rotation		// <-5/50,-10/10,-5/5>
transform {
	rotate <0,0,0>
	}

// Right leg rotations
define right_thigh_rotation		// <-90/15,0/-45,45/0>
transform {
	rotate <0,0,0>
	}
define right_calf_rotation		// <0/160,0,0>
transform {
	rotate <0,0,0>
	}
define right_foot_rotation		// <-5/50,-10/10,-5/5>
transform {
	rotate <0,0,0>
	}
