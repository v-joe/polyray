// This object is only included in the very first and very last frames,
// the first to generate depth information, and the last for the full
// image.  The sphere is not included when rendering the
if (frame < 7)
   define rsphere_tex matte_white
else
   define rsphere_tex
   texture {
      special surface {
	 color environment_map(reflect(I, N),
			       environment("out001.tga", "out002.tga",
					   "out003.tga", "out004.tga",
					   "out005.tga", "out006.tga"))
	 ambient 0.8
	 diffuse white, 0.2
	 specular white, 0.4
	 microfacet Phong 10
	 }
      }

// The reflective sphere itself.  We only need it in the scene when we are
// building the depth map for the light and when building the final image.
if (frame == 0 || frame == 7)
object {
   sphere sphere_location, 2
   rsphere_tex
   }
