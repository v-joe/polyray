// Woman, version 3 for Polyray, created with CSG.
// Created by Will Wagner.
// Thanks to Ken Boyce, Alexander Enzmann, and David Mason for their help.

include "rots.inc"

define woman
object {
	object {	// torso
		((object{sphere<0,2,0>,1.2 scale<1.665,1,1>}+
		object{cylinder<0,2,0>,<0,1.75,0>,1.2 scale<1.665,1,1>}+
		object{parabola<0,-1,0>,<0,1.75,0>,1.2 scale<1.665,1,1>})*
		object{box<1.6,3.5,-1.5>,<-1.6,-1,3>})+
		object{lathe 2,<0,0,1>,7,<0.8,0>,<0.725,0.2>,<0.65,0.5>,<0.55,0.7>,
			<0.35,0.9>,<0.15,0.975>,<0,1> shear 0,0,0,0,0,-0.65 rotate<5,10,0>
			translate<0.75,1.7,0.8> root_solver Sturm}+
		object{lathe 2,<0,0,1>,7,<0.8,0>,<0.725,0.2>,<0.65,0.5>,<0.55,0.7>,
			<0.35,0.9>,<0.15,0.975>,<0,1> shear 0,0,0,0,0,-0.65 rotate<5,-10,0>
			translate<-0.75,1.7,0.8> root_solver Sturm}+
		object {	//neck
			object{sphere<0,0,0>,0.575}+
			object{cylinder<0,0,0>,<0,1.5,0.45>,0.5}+
			object{sphere<0,1.5,0.45>,0.675}+
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
			object{sphere<0,0,0>,0.45}+
			object{cone<0,0,0>,0.45,<0,-2.4,0>,0.4}+
			object {	// right lower arm
				object{sphere<0,0,0>,0.4}+
				object{cone<0,0,0>,0.35,<0,-2.15>,0.2}+
				object{sphere<0,-2.15,0>,0.2}+
				object {	// right hand
					object{cone<0,0,0>,0.075,<0,-0.25,0>,0.15 scale<0.85,1,2.5>}+
					object{cylinder<0,-0.25,0>,<0,-0.8,0>,0.15 scale<0.85,1,2.5>}+
					object{sphere<0,-0.8,0>,0.15 scale<0.85,1,2.5>}+
					object {	// right fingers
						cone<0,0,0>,0.15,<0,-0.65,0>,0.1
						scale <1,1,3>
						right_finger_rotation
						translate<0,-0.8,0>
						}
					right_hand_rotation
					translate<0,-2.15,0>
					}
				right_lower_rotation
				translate<0,-2.4,0>
				}
			scale<0.925,1,1>
			right_upper_rotation
			translate<2.1,2.35,0>
			}+
		object {	// left upper arm
			object{sphere<0,0,0>,0.45}+
			object{cone<0,0,0>,0.45,<0,-2.4,0>,0.4}+
			object {	// left lower arm
				object{sphere<0,0,0>,0.4}+
				object{cone<0,0,0>,0.35,<0,-2.15>,0.2}+
				object{sphere<0,-2.15,0>,0.2}+
				object {	// left hand
					object{cone<0,0,0>,0.075,<0,-0.25,0>,0.15 scale<0.85,1,2.5>}+
					object{cylinder<0,-0.25,0>,<0,-0.8,0>,0.15 scale<0.85,1,2.5>}+
					object{sphere<0,-0.8,0>,0.15 scale<0.85,1,2.5>}+
					object {	// left fingers
						cone<0,0,0>,0.15,<0,-0.65,0>,0.1
						scale <0.85,1,2.5>
						left_finger_rotation
						translate<0,-0.8,0>
						}
					left_hand_rotation
					translate<0,-2.15,0>
					}
				left_lower_rotation
				translate<0,-2.4,0>
				}
			scale<0.925,1,1>
			left_upper_rotation
			translate<-2.1,2.35,0>
			}
		torso_rotation
		}+
	object {	// abdomen
		object{lathe 2,<0,1,0>,6,<0,1>,<-0.5,1>,<-1,0.75>,<-1,0>,<-0.25,-1>,<0,-1>
			root_solver Sturm scale<1.75,1,1>}+
		object {	// right thigh
			object{sphere<0,0,0>,0.8625}+
			object{cone<0,0,0>,0.85,<0,-2.75,0>,0.7}+
			object {	// right calf
				object{sphere<0,0,0>,0.7}+
				object{cone<0,0,0>,0.675,<0,-2.6,0>,0.4}+
				object {	// right foot
					object{disc<0,-0.4,0.875>,<0,0.9095,0.4158>,1.15}*
					object{disc<0,-0.85,1.75>,<0,0,1>,0.6}*
					object{disc<0,-0.75,0.525>,<0,-0.9919,-0.1266>,1.4}*
					object{disc<0,-0.3,-0.3>,<0,0.7071,-0.7071>,0.55}*
					object{disc<0.325,-0.85,0.525>,<0.9339,0.3512,-0.0673>,1.35}*
					object{disc<-0.325,-0.85,0.525>,<-0.9339,0.3512,-0.0673>,1.35}*
					(object{sphere<0,-0.85,1.2>,0.55}+
					object{cone<0,-0.85,1.2>,0.55,<0,-0.6,-0.6>,1})+
					object{sphere<0,0,0>,0.4}+
					object{sphere<0,-0.6,-0.6>,0.3}
					right_foot_rotation
					translate<0,-2.6,0>
					}
				right_calf_rotation
				translate<0,-2.75,0>
				}
			scale<0.925,1,1>
			right_thigh_rotation
			translate<0.9,-0.25,0>
			}+
		object {	// left thigh
			object{sphere<0,0,0>,0.8625}+
			object{cone<0,0,0>,0.85,<0,-2.75,0>,0.7}+
			object {	// left calf
				object{sphere<0,0,0>,0.7}+
				object{cone<0,0,0>,0.675,<0,-2.6,0>,0.4}+
				object {	// left foot
					object{disc<0,-0.4,0.875>,<0,0.9095,0.4158>,1.15}*
					object{disc<0,-0.85,1.75>,<0,0,1>,0.6}*
					object{disc<0,-0.75,0.525>,<0,-0.9919,-0.1266>,1.4}*
					object{disc<0,-0.3,-0.3>,<0,0.7071,-0.7071>,0.55}*
					object{disc<0.325,-0.85,0.525>,<0.9339,0.3512,-0.0673>,1.35}*
					object{disc<-0.325,-0.85,0.525>,<-0.9339,0.3512,-0.0673>,1.35}*
					(object{sphere<0,-0.85,1.2>,0.55}+
					object{cone<0,-0.85,1.2>,0.55,<0,-0.6,-0.6>,1})+
					object{sphere<0,0,0>,0.4}+
					object{sphere<0,-0.6,-0.6>,0.3}
					left_foot_rotation
					translate<0,-2.6,0>
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
	}
