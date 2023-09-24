//
// Slight adaptation of the CSG man developed by Will Wagner
//
define person
object {
   // Upper body
   object {
      // Head and neck
      object {
	 // Head
	 object {
	    object { sphere <0,-0.15,0.5>, 1 } +
	    ( object { cylinder <0,-0.15,-1.5>, <0,-0.15,0.5>, 1 } *
	      object { sphere <0,-1.025,0.5>, 1.375 } )
	    rotate headrot
	    translate <0,-0.45,1.5>
	    } +

	 // Neck
	 object { sphere <0,0,0>, 0.575 } +
	 object { cylinder <0,0,0>, <0,-0.45,1.5>, 0.575 } +
	 object { sphere <0,-0.45,1.5>, 0.675 }
	 rotate neckrot
	 translate <0,0.3,2.6>
	 } +

      // Right arm
      object {
	 // Lower arm
	 object {
	    // Hand
	    object {
	       // Fingers
	       object {
		  object {
		     cone <0,0,0>, 0.15, <0,0,-0.65>, 0.1
		     scale <1,3.5,1>
		     }
		  rotate rfingerrot
		  translate <0,0,-0.65>
		  } +

	       // Wrist & metacarpals
	       object { cylinder <0.15,0,0>, <-0.15,0,0>, 0.65 scale <1,1,0.5>}+
	       object { disc <0.15,0,0>, <1,0,0>, 0.65 scale <1,1,0.5> } +
	       object { disc <-0.15,0,0>, <-1,0,0>, 0.65 scale <1,1,0.5> } +
	       object { cylinder <0,0,0>, <0,0,-0.65>, 0.15 scale <1,3.5,1> } +
	       object { sphere <0,0,-0.65>, 0.15 scale<1,3.5,1> }
	       rotate rhandrot
	       translate <0,0,-2.25>
	       } +

	    // Elbow & humerus
	    object { sphere <0,0,0>, 0.6 } +
	    object { cone <0,0,0>, 0.5, <0,0,-2.25>, 0.2 }
	    rotate ralowerrot
	    translate <0,0,-2.575>
	    } +

	 // Shoulder & upper arm
	 object { sphere <0,0,0>, 0.75 } +
	 object { cone <0,0,0>, 0.65, <0,0,-2.575>, 0.5 }
	 rotate raupperrot
	 translate <2.65,0,2.25>
	 } +

      // Left arm
      object {
	 // Lower arm
	 object {
	    // Hand
	    object {
	       // Fingers
	       object {
		  object {
		     cone <0,0,0>, 0.15, <0,0,-0.65>, 0.1
		     scale<1,3.5,1>
		     }
		  rotate lfingerrot
		  translate <0,0,-0.65>
		  } +

	       // Wrist & metacarpals
	       object { cylinder <0.15,0,0>, <-0.15,0,0>, 0.65 scale <1,1,0.5>}+
	       object { disc <0.15,0,0>, <1,0,0>, 0.65 scale <1,1,0.5> } +
	       object { disc <-0.15,0,0>, <-1,0,0>, 0.65 scale <1,1,0.5> } +
	       object { cylinder <0,0,0>, <0,0,-0.65>, 0.15 scale <1,3.5,1> } +
	       object { sphere <0,0,-0.65>, 0.15 scale <1,3.5,1> }
	       rotate lhandrot
	       translate <0,0,-2.25>
	       } +

	    // Elbow & humerus
	    object { sphere <0,0,0>, 0.6 } +
	    object { cone <0,0,0>, 0.5, <0,0,-2.25>, 0.2 }
	    rotate lalowerrot
	    translate <0,0,-2.575>
	    } +

	 // Shoulder & upper arm
	 object { sphere <0,0,0>, 0.75 } +
	 object { cone <0,0,0>, 0.65, <0,0,-2.575>, 0.5 }
	 rotate laupperrot
	 translate <-2.65,0,2.25>
	 } +

      // Torso
      object { sphere <0,0,2>, 1.2 scale <1.85,1,1> } +
      object { cylinder <0,0,2>, <0,0,1.75>, 1.2 scale <1.85,1,1> } +
      object { parabola <0,0,-1>, <0,0,1.75>, 1.2 scale <1.85,1,1> }
      rotate torsorot
      } +

   // Lower body
   object {
      // Right leg
      object {
	 // Lower leg
	 object {
	    // Foot
	    object {
	       object { sphere <0,0,0>, 0.5 } +
	       object { cylinder <0,-0.25,-0.5>, <0,-0.25,-0.85>, 0.5
			scale <1,3.5,1> } +
	       ( object { sphere <0,-0.25,-0.5>, 0.5 scale <1,3.5,1> } *
		 object { disc <0,-0.25,-0.5>, <0,0,1>, 5 } )
	       rotate rfootrot
	       translate <0,0,-2.75>
	       } +

	    // Knee & Tibia
	    object { sphere <0,0,0>, 0.85 } +
	    object { cone <0,0,0>, 0.75, <0,0,-2.75>, 0.4 }
	    rotate rcalfrot
	    translate <0,0,-2.75>
	    } +
	 
	 // Hip joint
	 object { sphere <0,0,0>, 0.85 } +

	 // Femur & quads
	 object { cylinder <0,0,0>, <0,0,-2.75>, 0.8 }
	 rotate rthighrot
	 translate<1,0,-0.25>
	 } +

      // Left leg
      object {
	 // Lower leg
	 object {
	    // Foot
	    object {
	       object { sphere <0,0,0>, 0.5 } +
	       object { cylinder <0,-0.25,-0.5>, <0,-0.25,-0.85>, 0.5
			scale <1,3.5,1> } +
	       ( object { sphere <0,-0.25,-0.5>, 0.5 scale <1,3.5,1> } *
		 object { disc <0,-0.25,-0.5>, <0,0,1>, 5 } )
	       rotate lfootrot
	       translate <0,0,-2.75>
	       } +

	    // Knee & Tibia
	    object { sphere <0,0,0>, 0.85 } +
	    object { cone <0,0,0>, 0.75, <0,0,-2.75>, 0.4 }
	    rotate lcalfrot
	    translate <0,0,-2.75>
	    } +

	 // Hip joint
	 object { sphere <0,0,0>, 0.85 } +

	 // Femur & quads
	 object { cylinder <0,0,0>, <0,0,-2.75>, 0.8 }
	 rotate lthighrot
	 translate<-1,0,-0.25>
	 } +

      // Abdomen
      object {
	 lathe 2, <0,1,0>, 6,
	    <0,1,0>, <0.5,1,0>, <1,0.75,0>, <1,0,0>, <0.25,-1,0>, <0,-1,0>
	 root_solver Sturm
	 rotate<90,0,0>
	 scale<1.85,1,1>
	 }
      rotate abdomenrot
      translate<0,0,-1.75>
      }
   rotate <-90, 180, 0>
   translate <0, 8.35, 0>
   }
