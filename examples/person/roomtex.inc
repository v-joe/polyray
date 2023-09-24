define position_objectx     1
define position_cylindrical 3
define lookup_sawtooth 1

define light_wood <0.6, 0.3, 0.2>
define median_wood <0.4, 0.12, 0.06>
define dark_wood <0.2, 0.06, 0.005>
define reddish_wood
texture {
   noise surface {
      ambient 0.2
      diffuse 0.7
      specular white, 0.5
      microfacet Reitz 10
      reflection white, 0.15
      position_fn position_cylindrical
      position_scale 1
      lookup_fn lookup_sawtooth
      octaves 2
      turbulence 1
      color_map(
         [0.0, 0.8, light_wood, light_wood]
         [0.8, 1.0, median_wood, dark_wood])
      }
   scale <0.3, 0.3, 0.3>
   }

// The standard sort of marble texture
define green_marble
texture {
   noise surface {
      ambient 0.3
      diffuse 0.8
      specular 0.3
      reflection white, 0.3
      microfacet Reitz 5
      position_fn position_objectx
      lookup_fn lookup_sawtooth
      octaves 3
      turbulence 3
      color_map(
         [0.0, 0.8, aquamarine, 0.6 * aquamarine]
         [0.8, 1.0, 0.6 * aquamarine, <0.1, 0.1, 0.1>])
      }
   scale <0.5, 0.5, 0.5>
   }

// Nice blue agate texture
define sapphire_agate
texture {
   noise surface {
      ambient 0.5
      diffuse 0.7

      position_fn position_objectx
      position_scale 1.1
      lookup_fn lookup_sawtooth
      octaves 3
      turbulence 2

      color_map(
         [0.0, 0.3, <0, 0, 0.9>, <0, 0, 0.8>]
         [0.3, 1,   <0, 0, 0.8>, <0, 0, 0.4>])
      }
   scale <10, 10, 10>
   }

define orange_agate
texture {
   noise surface {
      ambient 0.5
      diffuse 0.7

      position_fn position_objectx
      position_scale 1.1
      lookup_fn lookup_sawtooth
      octaves 3
      turbulence 2

      color_map(
         [0.0, 0.3, <1, 0.2, 0.3>, 0.8 * <1, 0.2, 0.3>]
         [0.3, 1,   <1, 0.2, 0.3>, 0.2 * <1, 0.2, 0.3>])
      }
   scale <10, 10, 10>
   }

define waxed_white
texture {
   surface {
      color white
      ambient 0.1
      diffuse 0.6
      specular white, 0.3
      reflection white, 0.3
      microfacet Phong 10
      }
   }
