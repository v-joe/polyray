//
// Cratered landscape.
//
// A bump function was added to the texture to add to the apparent level of
// detail.  One way to reduce antialiasing in the render of this object is 
// to remove the bump.  This will give a flat shaded surface at distance.
//
// This object is composed of two objects: a flat plane (made from a disc)
// and a height field.  The plane is clipped by a box so that the height
// field can be placed into the middle.  The height field is clipped by
// a cylinder in the center so that the UFO launch pad can be inserted with
// out odd pieces of the hf showing up.
//
// The height field was generated using a paint program to make circles
// and lines.  An image processing program was then used to "blur" the
// starting image.  This gave smoother edges to the crater and rift edges.
//
static define craters
object {
  object {
	object { disc <0, 0, 0>, <0, 1, 0>, 1.5, 500 }
      & ~object { box <-128, -100, -128>, <128, 100, 128> }
      } +
   object {
      object {
	 smooth_height_field "crater2.jpg"
	 translate <-0.5, -48, -0.5>
	 rotate <0, 180, 0>
	 scale <256, 0.04, 256>
	 }
      & ~object { cylinder <0, -100, 0>, <0, 100, 0>, 1.5 }
      }
   /*
   texture {
      special matte {
	 color white
	 normal N + (dnoise(P, 3) - white/2)
	 }
      }
   */
   matte_white
   translate <0, -1, 0>
   }
