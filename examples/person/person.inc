// A person created with CSG.
// Created by Will Wagner.
// Thanks to Alexander Enzmann and Ken Boyce for their help.

include "rots.inc"

define head				//origin: top of neck
object {
	object{sphere<0,-0.15,0.5>,1}+
	(object{cylinder<0,-0.15,-1.5>,<0,-0.15,0.5>,1}*
	object{sphere<0,-1.025,0.5>,1.375})
	head_rotation
	translate<0,-0.45,1.5>
	}
define neck				//origin: base of neck
object {
	head+
	object{sphere<0,0,0>,0.575}+
	object{cylinder<0,0,0>,<0,-0.45,1.5>,0.575}+
	object{sphere<0,-0.45,1.5>,0.675}
	scale <0.85,1,1>
	neck_rotation
	translate<0,0.3,2.6>
	}
define r_fingers			//origin: palm/knuckles
object {
	object{cone<0,0,0>,0.15,<0,0,-0.65>,0.1 scale<1,3,1>}
	right_finger_rotation
	translate <0,0,-0.8>
	}
define l_fingers			//origin: palm/knuckles
object {
	object{cone<0,0,0>,0.15,<0,0,-0.65>,0.1 scale<1,3,1>}
	left_finger_rotation
	translate <0,0,-0.8>
	}
define r_hand				//origin: wrist
object {
	r_fingers+
	object{cone<0,0,0>,0.075,<0,0,-0.25>,0.15 scale<1,3,1>}+
	object{cylinder<0,0,-0.25>,<0,0,-0.8>,0.15 scale<1,3,1>}+
	object{sphere<0,0,-0.8>,0.15 scale<1,3,1>}
	right_hand_rotation
	translate<0,0,-2.25>
	}
define l_hand				//origin: wrist
object {
	l_fingers+
	object{cone<0,0,0>,0.075,<0,0,-0.25>,0.15 scale<1,3,1>}+
	object{cylinder<0,0,-0.25>,<0,0,-0.8>,0.15 scale<1,3,1>}+
	object{sphere<0,0,-0.8>,0.15 scale<1,3,1>}
	left_hand_rotation
	translate<0,0,-2.25>
	}
define r_lower_arm			//origin: elbow
object {
	r_hand+
	object{sphere<0,0,0>,0.5}+
	object{cone<0,0,0>,0.45,<0,0,-2.25>,0.2}+
	object{sphere<0,0,-2.25>,0.2}
	right_lower_rotation
	translate<0,0,-2.575>
	}
define l_lower_arm			//origin: elbow
object {
	l_hand+
	object{sphere<0,0,0>,0.5}+
	object{cone<0,0,0>,0.45,<0,0,-2.25>,0.2}+
	object{sphere<0,0,-2.25>,0.2}
	left_lower_rotation
	translate<0,0,-2.575>
	}
define r_upper_arm			//origin: shoulder
object {
	r_lower_arm+
	object{sphere<0,0,0>,0.55}+
	object{cone<0,0,0>,0.55,<0,0,-2.575>,0.5}
	scale <0.925,1,1>
	right_upper_rotation
	translate<2.35,0,2.25>
	}
define l_upper_arm			//origin: shoulder
object {
	l_lower_arm+
	object{sphere<0,0,0>,0.55}+
	object{cone<0,0,0>,0.55,<0,0,-2.575>,0.5}
	scale <0.925,1,1>
	left_upper_rotation
	translate<-2.35,0,2.25>
	}
define torso				//origin: lower back
object {
	neck+
	r_upper_arm+
	l_upper_arm+
	object{sphere<0,0,2>,1.2 scale<1.665,1,1>}+
	object{cylinder<0,0,2>,<0,0,1.75>,1.2 scale<1.665,1,1>}+
	object{parabola<0,0,-1>,<0,0,1.75>,1.2 scale<1.665,1,1>}
	torso_rotation
	}
define r_foot				//origin: ankle
object {
	object{sphere<0,0,0>,0.4}+
	object{sphere<0,0.45,-0.5>,0.4}+
	object{polygon 4,<0.35,0.6,-0.6>,<-0.35,0.6,-0.6>,<-0.6,-2.25,-0.9>,
		<0.6,-2.25,-0.9>}+
	object{polygon 4,<-0.6,-2.25,-0.8>,<0.6,-2.25,-0.8>,<0.1,0,0>,
		<-0.1,0,0>}+
	object{polygon 4,<-0.35,0.6,-0.6>,<-0.6,-2.25,-0.9>,<-0.6,-2.25,-0.8>,
		<-0.1,0,0>}+
	object{polygon 4,<0.35,0.6,-0.6>,<0.6,-2.25,-0.9>,<0.6,-2.25,-0.8>,
		<0.1,0,0>}+
	object{polygon 4,<0.1,0,0>,<-0.1,0,0>,<-0.35,0.6,-0.6>,<0.35,0.6,-0.6>}+
	object{polygon 4,<-0.6,-2.25,-0.9>,<-0.6,-2.25,-0.8>,<0.6,-2.25,-0.8>,
		<0.6,-2.25,-0.9>}
	right_foot_rotation
	translate<0,0,-2.75>
	}
define l_foot				//origin: ankle
object {
	object{sphere<0,0,0>,0.4}+
	object{sphere<0,0.45,-0.5>,0.4}+
	object{polygon 4,<0.35,0.6,-0.6>,<-0.35,0.6,-0.6>,<-0.6,-2.25,-0.9>,
		<0.6,-2.25,-0.9>}+
	object{polygon 4,<-0.6,-2.25,-0.8>,<0.6,-2.25,-0.8>,<0.1,0,0>,
		<-0.1,0,0>}+
	object{polygon 4,<-0.35,0.6,-0.6>,<-0.6,-2.25,-0.9>,<-0.6,-2.25,-0.8>,
		<-0.1,0,0>}+
	object{polygon 4,<0.35,0.6,-0.6>,<0.6,-2.25,-0.9>,<0.6,-2.25,-0.8>,
		<0.1,0,0>}+
	object{polygon 4,<0.1,0,0>,<-0.1,0,0>,<-0.35,0.6,-0.6>,<0.35,0.6,-0.6>}+
	object{polygon 4,<-0.6,-2.25,-0.9>,<-0.6,-2.25,-0.8>,<0.6,-2.25,-0.8>,
		<0.6,-2.25,-0.9>}
	left_foot_rotation
	translate<0,0,-2.75>
	}
define r_calf				//origin: knee
object {
	r_foot+
	object{sphere<0,0,0>,0.8}+
	object{cone<0,0,0>,0.75,<0,0,-2.75>,0.4}
	right_calf_rotation
	translate<0,0,-2.75>
	}
define l_calf				//origin: knee
object {
	l_foot+
	object{sphere<0,0,0>,0.8}+
	object{cone<0,0,0>,0.75,<0,0,-2.75>,0.4}
	left_calf_rotation
	translate<0,0,-2.75>
	}
define r_thigh				//origin: hip
object {
	r_calf+
	object{sphere<0,0,0>,0.825}+
	object{cylinder<0,0,0>,<0,0,-2.75>,0.8}
	scale <0.925,1,1>
	right_thigh_rotation
	translate<0.9,0,-0.25>
	}
define l_thigh				//origin: hip
object {
	l_calf+
	object{sphere<0,0,0>,0.825}+
	object{cylinder<0,0,0>,<0,0,-2.75>,0.8}
	scale <0.925,1,1>
	left_thigh_rotation
	translate<-0.9,0,-0.25>
	}
define abdomen				//origin: center of hips
object {
	r_thigh+
	l_thigh+
	object {lathe 2,<0,1,0>,6,<0,1,0>,<0.5,1,0>,<1,0.75,0>,<1,0,0>,
		<0.25,-1,0>,<0,-1,0> rotate<90,0,0> root_solver Sturm
		scale<1.665,1,1>}
	abdomen_rotation
	translate<0,0,-1.75>
	}
define person
object {
	torso+
	abdomen
	}