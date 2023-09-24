// Person, version 3 for Polyray, created with CSG.
// Created by Will Wagner.
// Thanks to Ken Boyce and Alexander Enzmann for their help.
// Test to see if combined form is faster

include "rots.inc"

define man
object {
	object {	// torso
		object{sphere<0,2,0>,1.2 scale<1.665,1,1>}+
		object{cylinder<0,2,0>,<0,1.75,0>,1.2 scale<1.665,1,1>}+
		object{parabola<0,-1,0>,<0,1.75,0>,1.2 scale<1.665,1,1>}+
		object {	// neck
			object{cone<0,0,0>,0.5,<0,1.5,0.45>,0.6}+
			object{sphere<0,1.5,0.45>,0.65}+
			object {	// head
				object{sphere<0,0.5,0.15>,1}+
				(object{cylinder<0,-1.5,0.15>,<0,0.5,0.15>,1}*
				object{sphere<0,0.5,1.025>,1.375})
				head_rotation
				translate<0,1.5,0.45>
				}
			scale<0.85,1,1>
			neck_rotation
			translate<0,2.6,-0.3>
			}+
		object {	// right upper arm
			object{sphere<0,0,0>,0.55}+
			object{cone<0,0,0>,0.55,<0,-2.575,0>,0.5}+
			object {	// right lower arm
				object{sphere<0,0,0>,0.5}+
				object{cone<0,0,0>,0.45,<0,-2.25>,0.2}+
				object{sphere<0,-2.25,0>,0.2}+
				object {	// right hand
					object{cone<0,0,0>,0.075,<0,-0.25,0>,0.15 scale<1,1,3>}+
					object{cylinder<0,-0.25,0>,<0,-0.8,0>,0.15 scale<1,1,3>}+
					object{sphere<0,-0.8,0>,0.15 scale<1,1,3>}+
					object {	// right fingers
						cone<0,0,0>,0.15,<0,-0.65,0>,0.1
						scale <1,1,3>
						right_finger_rotation
						translate<0,-0.8,0>
						}
					right_hand_rotation
					translate<0,-2.25,0>
					}
				right_lower_rotation
				translate<0,-2.575,0>
				}
			scale<0.925,1,1>
			right_upper_rotation
			translate<2.35,2.25,0>
			}+
		object {	// left upper arm
			object{sphere<0,0,0>,0.55}+
			object{cone<0,0,0>,0.55,<0,-2.575,0>,0.5}+
			object {	// left lower arm
				object{sphere<0,0,0>,0.5}+
				object{cone<0,0,0>,0.45,<0,-2.25>,0.2}+
				object{sphere<0,-2.25,0>,0.2}+
				object {	// left hand
					object{cone<0,0,0>,0.075,<0,-0.25,0>,0.15 scale<1,1,3>}+
					object{cylinder<0,-0.25,0>,<0,-0.8,0>,0.15 scale<1,1,3>}+
					object{sphere<0,-0.8,0>,0.15 scale<1,1,3>}+
					object {	// left fingers
						cone<0,0,0>,0.15,<0,-0.65,0>,0.1
						scale <1,1,3>
						left_finger_rotation
						translate<0,-0.8,0>
						}
					left_hand_rotation
					translate<0,-2.25,0>
					}
				left_lower_rotation
				translate<0,-2.575,0>
				}
			scale<0.925,1,1>
			left_upper_rotation
			translate<-2.35,2.25,0>
			}
		torso_rotation
		}+
	object {
		object {lathe 2,<0,1,0>,6,<0,1>,<-0.5,1>,<-1,0.75>,<-1,0>,<-0.25,-1>,<0,-1>
			root_solver Sturm scale<1.665,1,1>}+
			object {	// right thigh
				object{sphere<0,0,0>,0.825}+
				object{cylinder<0,0,0>,<0,-2.75,0>,0.8}+
				object {	// right calf
					object{sphere<0,0,0>,0.8}+
					object{cone<0,0,0>,0.75,<0,-2.75,0>,0.4}+
					object {	// right foot
						object{disc<0,-0.4,1.125>,<0,0.9422,0.3350>,1.35}*
						object{disc<0,-0.85,2.25>,<0,0,1>,0.65}*
						object{disc<0,-0.75,0.825>,<0,-0.9945,-0.1047>,1.60}*
						object{disc<0,-0.3,-0.3>,<0,0.7071,-0.7071>,0.55}*
						object{disc<0.475,-0.75,0.825>,<0.8910,0.4472,-0.0782>,1.5}*
						object{disc<-0.475,-0.75,0.825>,<-0.8910,0.4472,-0.0782>,1.5}*
						(object{sphere<0,-0.85,1.6>,0.65}+
						object{cone<0,-0.85,1.6>,0.65,<0,-0.6,-0.6>,1})+
						object{sphere<0,0,0>,0.4}+
						object{sphere<0,-0.5,-0.75>,0.4}
						right_foot_rotation
						translate<0,-2.75,0>
						}
					right_calf_rotation
					translate<0,-2.75,0>
					}
				scale<0.925,1,1>
				right_thigh_rotation
				translate<0.9,-0.25,0>
				}+
			object {	// left thigh
				object{sphere<0,0,0>,0.825}+
				object{cylinder<0,0,0>,<0,-2.75,0>,0.8}+
				object {	// left calf
					object{sphere<0,0,0>,0.8}+
					object{cone<0,0,0>,0.75,<0,-2.75,0>,0.4}+
					object {	// left foot
						object{disc<0,-0.4,1.125>,<0,0.9422,0.3350>,1.35}*
						object{disc<0,-0.85,2.25>,<0,0,1>,0.65}*
						object{disc<0,-0.75,0.825>,<0,-0.9945,-0.1047>,1.60}*
						object{disc<0,-0.3,-0.3>,<0,0.7071,-0.7071>,0.55}*
						object{disc<0.475,-0.75,0.825>,<0.8910,0.4472,-0.0782>,1.5}*
						object{disc<-0.475,-0.75,0.825>,<-0.8910,0.4472,-0.0782>,1.5}*
						(object{sphere<0,-0.85,1.6>,0.65}+
						object{cone<0,-0.85,1.6>,0.65,<0,-0.6,-0.6>,1})+
						object{sphere<0,0,0>,0.4}+
						object{sphere<0,-0.5,-0.75>,0.4}
						left_foot_rotation
						translate<0,-2.75,0>
						}
					left_calf_rotation
					translate<0,-2.75,0>
					}
				scale<0.925,1,1>
				left_thigh_rotation
				translate<-0.9,-0.25,0>
				}
			abdomen_rotation
			translate<0,-1.75,0>
			}
	personcolor
	}