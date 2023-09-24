// Definitions for a single UFO.  File by Alexander Enzmann
define UFO_Body_Size <5, 5, 5>

define UFO_Bounding_Scale <51, 41, 51>

// Define the length of the barrel of the guns.
define Gun_Length 10

define Gun_Assembly_Location1 <0, -20, 0>
define Gun_Assembly_Location2 <0,  20, 0>

// Define the locations of the two turrets on the gun.
define Gun_Location1 <-10, 0, 0>
define Gun_Location2 < 10, 0, 0>

// Define the location of the end of the barrel
define Barrel_Tip <0, 10, 0>

// Vary the length and width of the blast flame for a good blasting effect.
define Blast_Length 10
define Blast_Width 2.5

// Define the size of the chamber of the gun.
define Gun_Body_Size 3

// Define the size of the rod connecting the two barrels
define Gun_Rod_Connector_Size 2
define Transverse_Rod_Length 14
define Mounting_Rod_Length 20

// Define the size of the blast flame coming out of the barrel
define Gun_Blast_Scale <2, 0.4, 0.4>
define Gun_Blast_Offset <0, 11, 0>

// Very useful shape - used in this case to model the blast flame right
//  as a shot is fired.
define Piriform
object {
   polynomial 4*(x^4 - x^3) + y^2 + z^2
   bounding_box <0, -1, -1>, <1, 1, 1>
   shading_flags 32
   }

// Define a finite cylinder - this is useful for making rods of various
// lengths and widths using simple scaling.  This rod has length 1 and
// width 1.  The axis of the rod is along the Y-Axis.
define Rod
object {
   object { cylinder <0, 0, 0>, <0, 1, 0>, 1 }
   + object { disc <0, 0, 0>, <0, -1, 0>, 1 }
   + object { disc <0, 1, 0>, <0, 1, 0>, 1 }
   }

define Ground_Texture
texture {
   special matte {
      color color_map([0.0, 0.8, <0.8, 0.4, 0.2>, <1, 0.2, 0.4>]
		      [0.8, 1.0, <1, 0.2, 0.4>, <0.8, 0.4, 0.2>])[marble_fn]
      }
   scale <50, 50, 50>
   }

define Gun_Blast_Texture
texture { surface { color <1, 0.7, 0.3> ambient 1 diffuse 0 } }

// Define a gun barrel
define Gun_Barrel
object {
   object { cylinder <0, 0, 0>, <0, 10, 0>, 1 }
   + object { cylinder <0, 0, 0>, <0, 10, 0>, 0.6 }
   + object { disc <0, 0, 0>, <0, -1, 0>, 0.6, 1 }
   + object { disc <0, 10, 0>, <0, 1, 0>, 0.6, 1 }
   }

define Gun_Body
object { sphere <0, 0, 0>, Gun_Body_Size }

// Define the color of the gun blast
define Gun_Blast
Piriform {
   rotate <0, 0, -90>
   scale <Blast_Width, Blast_Length, Blast_Width>
   translate <0, Blast_Length, 0>
   Gun_Blast_Texture
   translate Barrel_Tip
   }

// The connecting frame for the guns is composed of two cylinders at
//  right angles, joined by a sphere.
define Gun_Connecting_Frame
object {
   // Transverse rod connecting the bodies of the two barrels }
   Rod {
      translate <0, -0.5, 0>
      scale <1, Transverse_Rod_Length, 1>
      rotate <0, 0, 90>
      }
    // Central connector from transverse rod to mounting rod }
    + object { sphere <0, 0, 0>, Gun_Rod_Connector_Size }
    Gun_Texture
    }

define Gun
object {
     Gun_Body { Gun_Texture }
   + Gun_Barrel { Gun_Texture }
   // + Gun_Blast { Gun_Blast_Texture }
   }

// Define the entire gun pod.  Two barrels, a cylinder connecting them,
// a sphere at the center of the cylinder for rotary control, and a cylinder
// extending up for a control rod.
define Gun_Barrels
object {
   // Place the two barrels of the gun }
     Gun { translate Gun_Location1 }
   + Gun { translate Gun_Location2 }
   // Transverse rod connecting the bodies of the two barrels
   + Rod {
      translate <0, -0.5, 0>
      scale <1, Transverse_Rod_Length, 1>
      rotate <0, 0, 90>
      Gun_Texture
      }
   // Central connector from transverse rod to mounting rod
   + object { sphere <0, 0, 0>, Gun_Rod_Connector_Size Gun_Texture }
   }

// Make UFO
define UFO_Body
object {
   // The main body of the UFO is made of scaled spheres.
     object { sphere <0, 0, 0>, 1 scale <10, 1, 10>
	      scale UFO_Body_Size Body_Tex1 }
   + object { sphere <0, 0, 0>, 1 scale <5, 2, 5>
	      scale UFO_Body_Size Body_Tex2 }
   + object { sphere <0, 0, 0>, 1 scale <3, 3, 3>
	      scale UFO_Body_Size Gun_Texture }
   // Mounting rod for the gun barrels.
   + Rod { scale <1, Mounting_Rod_Length, 1> Gun_Texture }
   }
