define position_cylindrical 3
define position_objectx 1
define lookup_sawtooth 1

//
// Some metallic surfaces
//
define brass_color  <0.71, 0.65, 0.26>
define copper_color <0.72, 0.45, 0.20>
define bronze_color <0.55, 0.47, 0.14>
define silver_color <0.90, 0.91, 0.98>

define brass
texture {
   surface {
      color brass_color
      ambient 0.2
      diffuse 0.7
      specular brass_color, 0.75
      microfacet Phong 10
      reflection 0.75
      }
   }

define chrome
texture {
   surface {
      color silver_color
      ambient 0.2
      diffuse 0.7
      specular silver_color, 0.75
      microfacet Phong 5
      reflection 0.2
      }
   }

define light_wood
texture {
   noise surface {
      ambient 0.2
      diffuse 0.7
      specular white, 0.5
      microfacet Reitz 10
      position_fn position_cylindrical
      position_scale 1
      lookup_fn lookup_sawtooth
      octaves 3
      turbulence 1
      color_map(
         [0.0, 0.8, <1, 0.72, 0.25>, <1, 0.72, 0.25>]
         [0.8, 1.0, <0.5, 0.5, 0.07>, <0.4, 0.4, 0.035>])
      }
   scale <0.4, 0.4, 0.4>
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
