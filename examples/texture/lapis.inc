//
// Adapted from the Persistence of Vision Raytracer
// Lapis Lazuli texture by Chad Eby
//

define lapis_layer1
color_map([0.0, 0.1, blue, midnight_blue]
          [0.1, 0.2, blue, midnight_blue]
          [0.2, 0.3, blue, midnight_blue]
          [0.3, 0.4, blue, midnight_blue]
          [0.4, 0.5, blue, midnight_blue]
          [0.5, 0.75, midnight_blue, slate_blue]
          [0.75, 0.9, slate_blue, gold]
	  [0.9, 0.97, midnight_blue, midnight_blue]
          [0.97, 1.0, slate_blue, gold])

define lapis_layer2
color_map([0.0, 0.4, black, 1.0, black, 1.0]
	  [0.4, 0.5, black, 1.0, white, 0.6]
	  [0.5, 0.7, white, 0.6, white, 0.6]
	  [0.7, 0.8, white, 0.6, black, 1.0]
	  [0.8, 0.97, black, 1.0, white/2, 0.5]
	  [0.97, 0.98, <0.2,0.2,1.0>, 0.5, blue, 0]
	  [0.98, 0.99, white, 0.5, white/2, 0.5]
	  [0.99, 1.0, white/2, 0.5, blue, 0])

static define lapis_layer2_fn sawtooth(x + 2*noise(P, 4))
static define lapis_layer1_fn noise(15*P, 5)
define Lapis_Lazuli
texture {
   layered
      texture { special shiny { color lapis_layer2[lapis_layer2_fn] } },
      texture { special shiny { color lapis_layer1[lapis_layer1_fn] } }
   }
