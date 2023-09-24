
// For any particular frame, we pull in the last one as the image map
if (frame == 0)
   define wall_image 0
else if (frame == 1)
   define wall_image image("out000.tga")
else if (frame == 2)
   define wall_image image("out001.tga")
else if (frame == 3)
   define wall_image image("out002.tga")
else if (frame == 4)
   define wall_image image("out003.tga")
else if (frame == 5)
   define wall_image image("out004.tga")
else if (frame == 6)
   define wall_image image("out005.tga")
else if (frame == 7)
   define wall_image image("out006.tga")
else if (frame == 8)
   define wall_image image("out007.tga")
else if (frame == 9)
   define wall_image image("out008.tga")
else if (frame == 10)
   define wall_image image("out009.tga")

//
// "wire_exper" returns 1 if we are on the visible part of the surface,
// and returns 0 if we are on the see-through part of the surface.
//
define density 0.5
define freq 0.5
define wall_exper (noise(freq*P, <2, 0.7, 3>) > density ? 1 : 0)
define wall_color  <0.71, 0.65, 0.26>
define lep_map
   color_map([0.0, 0.2, yellow, orange]
             [0.2, 0.4, orange, wall_color]
	     [0.4, 0.6, wall_color, brown]
	     [0.6, 0.8, brown, tan]
	     [0.8, 1.0, tan, sienna])
if (complexity == 0)
define wall_texture matte_white
else
define wall_texture
texture {
   layered
      texture {
	 special surface {
	    ambient wall_exper * wall_color, 0.2
	    diffuse wall_exper * wall_color, 0.5
	    specular wall_exper * wall_color, 0.5
	    microfacet Cook 5
	    transmission white, 1 - wall_exper, 1.0
	    // reflection white, wall_exper
	    }
	 },
      texture {
	 special matte {
	    color lep_map[leopard_fn]
	    }
	 }
      }

define diffusion_obj_map
   color_map([0.0, 0.2, white, white]
	     [0.2, 0.6, white, orange]
	     [0.6, 1.0, orange, coral])

define diffusion_obj_texture
   texture { special shiny { color diffusion_obj_map[noisy_onion_fn] } }

define Grnt0a
texture {
   noise surface {
      color white
      position_fn 1
      lookup_fn 1
      octaves 6
      turbulence 5
      ambient 0.2
      diffuse 0.6
      specular 0.3
      microfacet Reitz 5
      color_map(
	 [0.000, 0.153, <0.729, 0.502, 0.451>, 0.306,
			<0.769, 0.686, 0.592>, 0.792]
	 [0.153, 0.398, <0.769, 0.686, 0.592>,  0.792,
			<0.843, 0.753, 0.718>,  0.396]
	 [0.398, 0.559, <0.843, 0.753, 0.718>,  0.396,
			<0.780, 0.667, 0.561>,  0.976]
	 [0.559, 0.729, <0.780, 0.667, 0.561>,  0.976,
			<0.741, 0.659, 0.576>,  0.820]
	 [0.729, 1.001,  <0.741, 0.659, 0.576>, 0.820,
			 <0.729, 0.502, 0.451>,  0.306])
      }
   }

//-----Alpha Creams With Cracks
define Grnt9a
texture {
   noise surface {
      color white
      position_fn 1
      lookup_fn 1
      octaves 6
      turbulence 5
      ambient 0.2
      diffuse 0.6
      specular 0.3
      microfacet Reitz 5
   color_map(
      [0.000, 0.216, <0.812, 0.812, 0.812>, 0.835,
      <0.745, 0.843, 0.835>, 0.847]
      [0.216, 0.241, <0.745, 0.843, 0.835>, 0.847,
      <0.404, 0.337, 0.337>, 0.463]
      [0.241, 0.267, <0.404, 0.337, 0.337>, 0.463,
      <0.773, 0.729, 0.745>, 0.622]
      [0.267, 0.759, <0.773, 0.729, 0.745>, 0.622,
      <0.914, 0.843, 0.725>, 0.651]
      [0.759, 0.784, <0.914, 0.843, 0.725>, 0.651,
      <0.153, 0.133, 0.208>, 0.437]
      [0.784, 0.810, <0.153, 0.133, 0.208>, 0.437,
      <0.812, 0.812, 0.812>, 0.835]
      [0.810, 1.001, <0.812, 0.812, 0.812>, 0.835,
      <0.812, 0.812, 0.812>, 0.835])
   }
}

//---------- Transparent Green Moss Colors
define Grnt23a
texture {
   noise surface {
      color white
      position_fn 1
      lookup_fn 1
      octaves 6
      turbulence 5
      ambient 0.2
      diffuse 0.6
      specular 0.3
      microfacet Reitz 5
   color_map(
      [0.000, 0.168, <0.824, 0.725, 0.584>, 0.600,
      <0.514, 0.584, 0.533>, 0.600]
      [0.168, 0.301, <0.514, 0.584, 0.533>, 0.600,
      <0.298, 0.376, 0.318>, 0.600]
      [0.301, 0.398, <0.298, 0.376, 0.318>, 0.600,
      <0.263, 0.337, 0.282>, 0.700]
      [0.398, 0.558, <0.263, 0.337, 0.282>, 0.700,
      <0.431, 0.506, 0.451>, 0.600]
      [0.558, 0.655, <0.431, 0.506, 0.451>, 0.600,
      <0.529, 0.631, 0.471>, 0.500]
      [0.655, 0.735, <0.529, 0.631, 0.471>, 0.500,
      <0.333, 0.376, 0.318>, 0.700]
      [0.735, 0.823, <0.333, 0.376, 0.318>, 0.700,
      <0.298, 0.376, 0.318>, 0.600]
      [0.823, 0.876, <0.298, 0.376, 0.318>, 0.600,
      <0.416, 0.376, 0.318>, 0.500]
      [0.876, 0.929, <0.416, 0.376, 0.318>, 0.500,
      <0.416, 0.376, 0.318>, 0.600]
      [0.929, 1.001, <0.416, 0.376, 0.318>, 0.600,
      <0.824, 0.725, 0.584>, 0.700])
   }
}

define Grnt24a
texture {
   noise surface {
      color white
      position_fn 1
      lookup_fn 1
      octaves 6
      turbulence 5
      ambient 0.2
      diffuse 0.6
      specular 0.3
      microfacet Reitz 5
   color_map(
      [0.000, 0.053, <0.784, 0.627, 0.522>, 0.500,
      <0.784, 0.627, 0.624>, 0.500]
      [0.053, 0.263, <0.784, 0.627, 0.624>, 0.500,
      <0.824, 0.557, 0.376>, 0.500]
      [0.263, 0.281, <0.824, 0.557, 0.376>, 0.500,
      <0.643, 0.380, 0.376>, 0.500]
      [0.281, 0.325, <0.643, 0.380, 0.376>, 0.500,
      <0.839, 0.722, 0.722>, 0.500]
      [0.325, 0.711, <0.839, 0.722, 0.722>, 0.500,
      <0.784, 0.627, 0.522>, 0.500]
      [0.711, 0.798, <0.784, 0.627, 0.522>, 0.500,
      <0.769, 0.380, 0.376>, 0.500]
      [0.798, 0.895, <0.769, 0.380, 0.376>, 0.500,
      <0.824, 0.557, 0.376>, 0.500]
      [0.895, 0.982, <0.824, 0.557, 0.376>, 0.500,
      <0.784, 0.627, 0.522>, 0.500]
      [0.982, 1.001, <0.784, 0.627, 0.522>, 0.500,
      <0.784, 0.627, 0.522>, 0.500])
   }
}

//--------------------Crack & OverTint /Red---------------------
define Crack1
texture {
   noise surface {
      color white
      position_fn 1
      lookup_fn 1
      octaves 6
      turbulence 5
      ambient 0.1
      diffuse 0.2
      specular 0.3
      microfacet Reitz 5
   color_map(
      [0.00, 0.04, black,     0.6, black,   1.0]
      [0.04, 0.97, firebrick, 0.8, dimgray, 0.9] 
      [0.97, 1.00, black,     0.9, black,   1.0])
   }
}       

define clear_yellow
texture {
   surface {
      ambient 0
      diffuse 0
      transmission yellow, 0.9, 1
      }
   }

//------------ Moss Greens & Tan 
define Grnt25
texture {
   noise surface {
      color white
      position_fn 1
      lookup_fn 1
      octaves 6
      turbulence 5
      ambient 0.2
      diffuse 0.6
      specular 0.3
      microfacet Reitz 5
      color_map(
	 [0.000, 0.168, <0.824, 0.725, 0.584>, 0.000,
			<0.514, 0.584, 0.533>, 0.000]
	 [0.168, 0.301, <0.514, 0.584, 0.533>, 0.000,
			<0.298, 0.376, 0.318>, 0.000]
	 [0.301, 0.398, <0.298, 0.376, 0.318>, 0.000,
			<0.263, 0.337, 0.282>, 0.000]
	 [0.398, 0.558, <0.263, 0.337, 0.282>, 0.000,
			<0.431, 0.506, 0.451>, 0.000]
	 [0.558, 0.655, <0.431, 0.506, 0.451>, 0.000,
			<0.529, 0.631, 0.471>, 0.000]
	 [0.655, 0.735, <0.529, 0.631, 0.471>, 0.000,
			<0.333, 0.376, 0.318>, 0.000]
	 [0.735, 0.823, <0.333, 0.376, 0.318>, 0.000,
			<0.298, 0.376, 0.318>, 0.000]
	 [0.823, 0.876, <0.298, 0.376, 0.318>, 0.000,
			<0.416, 0.376, 0.318>, 0.000]
	 [0.876, 0.929, <0.416, 0.376, 0.318>, 0.000,
			<0.416, 0.376, 0.318>, 0.000]
	 [0.929, 1.001, <0.416, 0.376, 0.318>, 0.000,
			<0.824, 0.725, 0.584>, 0.000])
   }
}

//-------- White Tan & thin Reds
define Grnt29
texture {
   noise surface {
      color white
      position_fn 1
      lookup_fn 1
      octaves 6
      turbulence 5
      ambient 0.2
      diffuse 0.6
      specular 0.3
      microfacet Reitz 5
   color_map(
      [0.000, 0.053, <0.784, 0.627, 0.522>, 0.000,
      <0.784, 0.627, 0.624>, 0.000]
      [0.053, 0.263, <0.784, 0.627, 0.624>, 0.000,
      <0.824, 0.557, 0.376>, 0.000]
      [0.263, 0.281, <0.824, 0.557, 0.376>, 0.000,
      <0.643, 0.380, 0.376>, 0.000]
      [0.281, 0.325, <0.643, 0.380, 0.376>, 0.000,
      <0.839, 0.722, 0.722>, 0.000]
      [0.325, 0.711, <0.839, 0.722, 0.722>, 0.000,
      <0.784, 0.627, 0.522>, 0.000]
      [0.711, 0.798, <0.784, 0.627, 0.522>, 0.000,
      <0.769, 0.380, 0.376>, 0.000]
      [0.798, 0.895, <0.769, 0.380, 0.376>, 0.000,
      <0.824, 0.557, 0.376>, 0.000]
      [0.895, 0.982, <0.824, 0.557, 0.376>, 0.000,
      <0.784, 0.627, 0.522>, 0.000]
      [0.982, 1.001, <0.784, 0.627, 0.522>, 0.000,
      <0.784, 0.627, 0.522>, 0.000])
   }
}

//------------- Gray Granite with white patches
define Stone8
texture {
   layered
      Grnt0a { scale <2, 3, 2> rotate <0, 0, -30> },
      Grnt9a { scale <5, 3, 4> rotate <0, 0, 40> },
      texture { shiny { color white } }
   }


//------------- Peach & Cream Marble with orange veining
define Stone23
texture {
   layered
      clear_yellow,
      Crack1 { scale <1, 2, 1.5> rotate <0, 0, 40> },
      Grnt24a { scale <2, 1, 2> rotate <40, 0, 30> },
      Grnt29 { scale <1, 1, 2> rotate <40, 0, 30> }
   }

//------------- Green & Tan Moss Marble
define Stone24
texture {
   layered
      Grnt23a { scale <2, 1, 2> rotate <40, 0, -30> },
      Grnt25 { scale <1, 1, 2> rotate <0, 0, 50> rotate <20, 0, 30> }
   }

define jade
texture {
   noise surface {
      color white
      ambient 0.3
      diffuse 0.8
      specular 0.3
      microfacet Reitz 5
      position_fn position_objectx
      lookup_fn   lookup_sawtooth
      octaves     3
      turbulence  3
      color_map(
         [0.0, 0.8, <0.1, 0.6, 0.1>, <0.0, 0.3, 0.0>]
         [0.8, 1.0, <0.1, 0.6, 0.1>, <0.0, 0.2, 0.0>])
      }
    rotate <80, 0, 0>
    scale <0.5, 0.5, 0.5>
   }
