//
//  This set of textures is based on the original stones textures developed
//  for the "Persistence of Vision Raytracer".  It has been included with
//  the Polyray data archives with permission of the original artist.
//
//		  Stone textures by Mike Miller  1992
//
//  You can simply use the textures Stone1 -> Stone24, or you can develop
//  other similar textures by using the "Grnt-a" textures layered over the
//  Grnt textures.
//
define noisy_layers (sin(2*pi*(x + 2*fnoise(P, 3)))+1)/2 // Like marble
define noisy_layers2 sawtooth(x + 2*fnoise(P, 3))

define granite_fn_04 fnoise(4 * (P + 0.8 * dnoise(P, 1)), 6)
define granite_fn_05 fnoise(4 * (P + 1.0 * dnoise(P, 1)), 6)
define granite_fn_06 fnoise(4 * (P + 1.2 * dnoise(P, 1)), 6)
define granite_fn_07 fnoise(4 * (P + 1.4 * dnoise(P, 1)), 6)
define granite_fn fnoise(4 * P, 6)

define shiny_white
texture {
   surface {
      ambient white, 0.1
      diffuse white, 0.3
      specular white, 0.4
      microfacet Phong 3
      }
   }

define shiny_black
texture {
   surface {
      ambient black, 0.0
      diffuse white, 0.1
      specular white, 0.8
      microfacet Phong 3
      }
   }

define clear_yellow
texture {
   surface {
      color yellow
      ambient 0.1
      diffuse 0.3
      transmission yellow, 0.8, 1
      specular white, 0.8
      microfacet Phong 3
      }
   }

//--------- Gray  Tan with Rose
define Grnt0_map
      color_map(
	 [0.000, 0.153, <0.729, 0.502, 0.451>, <0.769, 0.686, 0.592>]
	 [0.153, 0.398, <0.769, 0.686, 0.592>, <0.843, 0.753, 0.718>]
	 [0.398, 0.559, <0.843, 0.753, 0.718>, <0.780, 0.667, 0.561>]
	 [0.559, 0.729, <0.780, 0.667, 0.561>, <0.741, 0.659, 0.576>]
	 [0.729, 1.001, <0.741, 0.659, 0.576>, <0.729, 0.502, 0.451>])

//----- Creamy Whites with yellow & light gray
define Grnt1_map
      color_map(
	 [0.000, 0.212, <0.898, 0.898, 0.851>, <0.969, 0.980, 0.875>]
	 [0.212, 0.424, <0.969, 0.980, 0.875>, <0.859, 0.859, 0.859>]
	 [0.424, 0.627, <0.859, 0.859, 0.859>, <0.992, 0.922, 0.659>]
	 [0.627, 0.881, <0.992, 0.922, 0.659>, <0.937, 0.965, 0.902>]
	 [0.881, 1.000, <0.937, 0.965, 0.902>, <0.898, 0.898, 0.851>])

//------- Deep Cream with light rose, yellow orchid & tan
define Grnt2_map
      color_map(
	 [0.000, 0.178, <0.863, 0.757, 0.596>, <0.925, 0.792, 0.714>]
	 [0.178, 0.356, <0.925, 0.792, 0.714>, <0.871, 0.702, 0.659>]
	 [0.356, 0.525, <0.871, 0.702, 0.659>, <0.992, 0.922, 0.659>]
	 [0.525, 0.729, <0.992, 0.922, 0.659>, <0.902, 0.812, 0.714>]
	 [0.729, 1.001, <0.902, 0.812, 0.714>, <0.863, 0.757, 0.596>])

//------- Warm tans olive & light rose with cream
define Grnt3_map
      color_map(
	 [0.000, 0.178, <0.831, 0.631, 0.569>, <0.925, 0.831, 0.714>]
	 [0.178, 0.356, <0.925, 0.831, 0.714>, <0.871, 0.702, 0.659>]
	 [0.356, 0.525, <0.871, 0.702, 0.659>, <0.831, 0.631, 0.569>]
	 [0.525, 0.729, <0.831, 0.631, 0.569>, <0.937, 0.882, 0.820>]
	 [0.729, 1.001, <0.937, 0.882, 0.820>, <0.831, 0.631, 0.569>])

//-------- Orchid sand & mouve
define Grnt4_map
      color_map(
	 [0.000, 0.178, <0.804, 0.569, 0.494>, <0.816, 0.725, 0.537>]
	 [0.178, 0.356, <0.816, 0.725, 0.537>, <0.820, 0.580, 0.522>]
	 [0.356, 0.525, <0.820, 0.580, 0.522>, <0.882, 0.725, 0.537>]
	 [0.525, 0.729, <0.882, 0.725, 0.537>, <0.855, 0.729, 0.584>]
	 [0.729, 1.001, <0.855, 0.729, 0.584>, <0.804, 0.569, 0.494>])

//------- Medium Mauve Med.Rose & deep cream
define Grnt5_map
      color_map(
	 [0.000, 0.178, <0.804, 0.569, 0.494>, <0.855, 0.729, 0.584>]
	 [0.178, 0.356, <0.855, 0.729, 0.584>, <0.667, 0.502, 0.478>]
	 [0.356, 0.525, <0.667, 0.502, 0.478>, <0.859, 0.624, 0.545>]
	 [0.525, 0.729, <0.859, 0.624, 0.545>, <0.855, 0.729, 0.584>]
	 [0.729, 1.001, <0.855, 0.729, 0.584>, <0.804, 0.569, 0.494>])

//--------- Med. Orchid Olive & Dark Tan "mud pie"
define Grnt6_map
      color_map(
	 [0.000, 0.153, <0.545, 0.380, 0.345>, <0.588, 0.475, 0.333>]
	 [0.153, 0.398, <0.588, 0.475, 0.333>, <0.675, 0.478, 0.404>]
	 [0.398, 0.559, <0.675, 0.478, 0.404>, <0.757, 0.635, 0.522>]
	 [0.559, 0.729, <0.757, 0.635, 0.522>, <0.659, 0.549, 0.443>]
	 [0.729, 1.001, <0.659, 0.549, 0.443>, <0.545, 0.380, 0.345>])

//------- Dark Orchid Olive & Dark Putty
define Grnt7_map
      color_map(
	 [0.000, 0.119,  <0.439, 0.310, 0.282>, <0.463, 0.369, 0.259>]
	 [0.119, 0.322,  <0.463, 0.369, 0.259>, <0.541, 0.369, 0.298>]
	 [0.322, 0.449,  <0.541, 0.369, 0.298>, <0.573, 0.424, 0.286>]
	 [0.449, 0.729,  <0.573, 0.424, 0.286>, <0.494, 0.396, 0.306>]
	 [0.729, 1.001,  <0.494, 0.396, 0.306>, <0.439, 0.310, 0.282>])

//--------- Rose & Light cream Yellows
define Grnt8_map
      color_map(
	 [0.000, 0.179, <0.843, 0.655, 0.655>, <0.886, 0.769, 0.627>]
	 [0.179, 0.368, <0.886, 0.769, 0.627>, <0.906, 0.820, 0.714>]
	 [0.368, 0.538, <0.906, 0.820, 0.714>, <0.851, 0.671, 0.671>]
	 [0.538, 0.846, <0.851, 0.671, 0.671>, <0.890, 0.792, 0.675>]
	 [0.846, 0.983, <0.890, 0.792, 0.675>, <0.827, 0.612, 0.612>]
	 [0.983, 1.001, <0.827, 0.612, 0.612>, <0.843, 0.655, 0.655>])

//--------- Light Steely Grays
define Grnt9_map
      color_map(
	 [0.000, 0.154, <0.894, 0.886, 0.886>, <0.745, 0.745, 0.753>]
	 [0.154, 0.308, <0.745, 0.745, 0.753>, <0.902, 0.902, 0.859>]
	 [0.308, 0.444, <0.902, 0.902, 0.859>, <0.729, 0.706, 0.694>]
	 [0.444, 0.615, <0.729, 0.706, 0.694>, <0.588, 0.592, 0.635>]
	 [0.615, 0.803, <0.588, 0.592, 0.635>, <0.608, 0.616, 0.659>]
	 [0.803, 1.001, <0.608, 0.616, 0.659>, <0.894, 0.886, 0.886>])

//--------- Gray Creams & lavender tans
define Grnt10_map
      color_map(
	 [0.000, 0.154, <0.890, 0.690, 0.690>, <0.996, 0.835, 0.737>]
	 [0.154, 0.308, <0.996, 0.835, 0.737>, <0.745, 0.635, 0.651>]
	 [0.308, 0.444, <0.745, 0.635, 0.651>, <0.733, 0.596, 0.557>]
	 [0.444, 0.615, <0.733, 0.596, 0.557>, <0.996, 0.835, 0.737>]
	 [0.615, 0.803, <0.996, 0.835, 0.737>, <0.765, 0.616, 0.659>]
	 [0.803, 1.001, <0.765, 0.616, 0.659>, <0.890, 0.690, 0.690>])

//--------- Creams & Grays Kakhi
define Grnt11_map
      color_map(
	 [0.000, 0.154, <0.800, 0.651, 0.557>, <0.996, 0.835, 0.737>]
	 [0.154, 0.308, <0.996, 0.835, 0.737>, <0.800, 0.651, 0.557>]
	 [0.308, 0.444, <0.800, 0.651, 0.557>, <0.694, 0.624, 0.604>]
	 [0.444, 0.615, <0.694, 0.624, 0.604>, <0.800, 0.651, 0.557>]
	 [0.615, 0.812, <0.800, 0.651, 0.557>, <0.725, 0.655, 0.651>]
	 [0.812, 1.001, <0.725, 0.655, 0.651>, <0.800, 0.651, 0.557>])

//--------- Tan Cream & Red Rose 
define Grnt12_map
      color_map(
	 [0.000, 0.154, <0.996, 0.969, 0.800>, <0.996, 0.682, 0.604>]
	 [0.154, 0.308, <0.996, 0.682, 0.604>, <0.906, 0.820, 0.714>]
	 [0.308, 0.444, <0.906, 0.820, 0.714>, <0.816, 0.631, 0.537>]
	 [0.444, 0.615, <0.816, 0.631, 0.537>, <0.890, 0.792, 0.675>]
	 [0.615, 0.812, <0.890, 0.792, 0.675>, <0.973, 0.627, 0.627>]
	 [0.812, 1.001, <0.973, 0.627, 0.627>, <0.996, 0.969, 0.800>])

//--------- Cream Rose orange 
define Grnt13_map
      color_map(
	 [0.000, 0.154, <0.996, 0.824, 0.780>, <0.996, 0.698, 0.624>]
	 [0.154, 0.308, <0.996, 0.698, 0.624>, <0.906, 0.675, 0.553>]
	 [0.308, 0.444, <0.906, 0.675, 0.553>, <0.996, 0.682, 0.604>]
	 [0.444, 0.615, <0.996, 0.682, 0.604>, <0.996, 0.824, 0.780>]
	 [0.615, 0.812, <0.996, 0.824, 0.780>, <0.973, 0.627, 0.627>]
	 [0.812, 1.001, <0.973, 0.627, 0.627>, <0.996, 0.824, 0.780>])

//--------- Cream Rose & light moss & light Violet
define Grnt14_map
      color_map(
	 [0.000, 0.154, <0.690, 0.612, 0.569>, <0.737, 0.596, 0.522>]
	 [0.154, 0.368, <0.737, 0.596, 0.522>, <0.776, 0.702, 0.624>]
	 [0.368, 0.538, <0.776, 0.702, 0.624>, <0.796, 0.678, 0.643>]
	 [0.538, 0.846, <0.796, 0.678, 0.643>, <0.690, 0.612, 0.569>]
	 [0.846, 0.932, <0.690, 0.612, 0.569>, <0.773, 0.612, 0.569>]
	 [0.932, 1.001, <0.773, 0.612, 0.569>, <0.690, 0.612, 0.569>])

//--------- Black with subtle chroma
define Grnt15_map
      color_map(
	 [0.000, 0.104, <0.161, 0.133, 0.118>, <0.110, 0.082, 0.071>]
	 [0.104, 0.252, <0.110, 0.082, 0.071>, <0.161, 0.133, 0.118>]
	 [0.252, 0.383, <0.161, 0.133, 0.118>, <0.000, 0.000, 0.000>]
	 [0.383, 0.643, <0.000, 0.000, 0.000>, <0.161, 0.133, 0.118>]
	 [0.643, 0.783, <0.161, 0.133, 0.118>, <0.220, 0.149, 0.137>]
	 [0.783, 0.922, <0.220, 0.149, 0.137>, <0.000, 0.000, 0.000>]
	 [0.922, 0.983, <0.000, 0.000, 0.000>, <0.220, 0.149, 0.137>]
	 [0.983, 1.001, <0.220, 0.149, 0.137>, <0.161, 0.133, 0.118>])

//----- White Cream & Peach
define Grnt16_map
      color_map(
	 [0.000, 0.316, <0.910, 0.788, 0.788>, <0.922, 0.914, 0.871>]
	 [0.316, 0.453, <0.922, 0.914, 0.871>, <0.894, 0.867, 0.780>]
	 [0.453, 0.624, <0.894, 0.867, 0.780>, <0.784, 0.788, 0.788>]
	 [0.624, 0.726, <0.784, 0.788, 0.788>, <0.851, 0.812, 0.741>]
	 [0.726, 0.863, <0.851, 0.812, 0.741>, <0.647, 0.655, 0.655>]
	 [0.863, 1.001, <0.647, 0.655, 0.655>, <0.910, 0.788, 0.788>])

//----- Bug Juice & Green
define Grnt17_map
      color_map(
	 [0.000, 0.303, <0.000, 0.239, 0.000>, <0.333, 0.294, 0.000>]
	 [0.303, 0.588, <0.333, 0.294, 0.000>, <0.000, 0.239, 0.341>]
	 [0.588, 0.790, <0.000, 0.239, 0.341>, <0.000, 0.020, 0.000>]
	 [0.790, 1.001, <0.000, 0.020, 0.000>, <0.000, 0.239, 0.000>])

//------------ Rose & cream yellow
define Grnt18_map
      color_map(
	 [0.000, 0.202, <1.000, 0.718, 0.541>, <0.890, 0.651, 0.612>]
	 [0.202, 0.298, <0.890, 0.651, 0.612>, <1.000, 0.820, 0.675>]
	 [0.298, 0.377, <1.000, 0.820, 0.675>, <0.890, 0.643, 0.612>]
	 [0.377, 0.465, <0.890, 0.643, 0.612>, <0.937, 0.729, 0.561>]
	 [0.465, 0.544, <0.937, 0.729, 0.561>, <0.878, 0.604, 0.565>]
	 [0.544, 0.640, <0.878, 0.604, 0.565>, <0.984, 0.780, 0.655>]
	 [0.640, 0.860, <0.984, 0.780, 0.655>, <1.000, 0.863, 0.635>]
	 [0.860, 0.982, <1.000, 0.863, 0.635>, <1.000, 0.765, 0.620>]
	 [0.982, 1.001, <1.000, 0.765, 0.620>, <1.000, 0.718, 0.541>])

//--------- Gray Marble with White feather Viens
define Grnt19_map
      color_map(
	 [0.00, 0.35, White,   White]
	 [0.35, 0.40, White,   DimGray]
	 [0.40, 1.00, DimGray, DimGray])

//--------- White Marble with Gray feather Viens
define Grnt20_map
      color_map(
	 [0.0, 0.3, black, white]
	 [0.2, 0.4, white, white]
	 [0.4, 0.6, white, white]
	 [0.6, 1.0, white, white])

//-------- Declare Green Colors
define g1 <0.26, 0.41, 0.31> //---Light Gray Green
define g2 <0.27, 0.34, 0.26> //---Med Gray Green
define g3 <0.13, 0.29, 0.28> //---Med Gray Aqua
define g4 <0.03, 0.18, 0.08> //---Dark Green

//--------- Green Jade
define Grnt21_map
      color_map(
	 [0.0, 0.1, White,     SeaGreen ]
	 [0.1, 0.3, SeaGreen,  g2       ]
	 [0.3, 0.5, g2,        DarkGreen]
	 [0.5, 0.7, DarkGreen, g4       ]
	 [0.7, 0.8, g4,        DarkGreen]
	 [0.8, 1.0, DarkGreen, DarkGreen])

//--------- Clear with White feather Viens ----- This one does contain Alpha
define Grnt22_map
      color_map(
	 [0.00, 0.07, White, White]
	 [0.07, 0.20, White, DimGray]
	 [0.20, 0.30, DimGray, 0.0, white, 1.0]
	 [0.30, 0.70, white,   1.0, white, 1.0]
	 [0.70, 1.00, white,   1.0, DimGray, 0.0])

//---------- Light Tan to Mouve 
define Grnt23_map
      color_map(
	 [0.000, 0.178, <0.831, 0.631, 0.569>, <0.925, 0.831, 0.714>]
	 [0.178, 0.356, <0.925, 0.831, 0.714>, <0.871, 0.702, 0.659>]
	 [0.356, 0.525, <0.871, 0.702, 0.659>, <0.831, 0.631, 0.569>]
	 [0.525, 0.729, <0.831, 0.631, 0.569>, <0.937, 0.882, 0.820>]
	 [0.729, 1.001, <0.937, 0.882, 0.820>, <0.831, 0.631, 0.569>])

//--------- Light Grays
define Grnt24_map
      color_map(
	 [0.000, 0.154, <0.894, 0.886, 0.886>, <0.745, 0.745, 0.753>]
	 [0.154, 0.308, <0.745, 0.745, 0.753>, <0.902, 0.902, 0.859>]
	 [0.308, 0.444, <0.902, 0.902, 0.859>, <0.729, 0.706, 0.694>]
	 [0.444, 0.615, <0.729, 0.706, 0.694>, <0.588, 0.592, 0.635>]
	 [0.615, 0.803, <0.588, 0.592, 0.635>, <0.608, 0.616, 0.659>]
	 [0.803, 1.001, <0.608, 0.616, 0.659>, <0.894, 0.886, 0.886>])

//------------ Moss Greens & Tan 
define Grnt25_map
      color_map(
	 [0.000, 0.168, <0.824, 0.725, 0.584>, <0.514, 0.584, 0.533>]
	 [0.168, 0.301, <0.514, 0.584, 0.533>, <0.298, 0.376, 0.318>]
	 [0.301, 0.398, <0.298, 0.376, 0.318>, <0.263, 0.337, 0.282>]
	 [0.398, 0.558, <0.263, 0.337, 0.282>, <0.431, 0.506, 0.451>]
	 [0.558, 0.655, <0.431, 0.506, 0.451>, <0.529, 0.631, 0.471>]
	 [0.655, 0.735, <0.529, 0.631, 0.471>, <0.333, 0.376, 0.318>]
	 [0.735, 0.823, <0.333, 0.376, 0.318>, <0.298, 0.376, 0.318>]
	 [0.823, 0.876, <0.298, 0.376, 0.318>, <0.416, 0.376, 0.318>]
	 [0.876, 0.929, <0.416, 0.376, 0.318>, <0.416, 0.376, 0.318>]
	 [0.929, 1.001, <0.416, 0.376, 0.318>, <0.824, 0.725, 0.584>])

//---------- Salmon with thin Green Viens
define Grnt26_map
   color_map(
      [0.000, 0.241, <0.973, 0.973, 0.976>, <0.973, 0.973, 0.976>]
      [0.241, 0.284, <0.973, 0.973, 0.976>, <0.600, 0.741, 0.608>]
      [0.284, 0.336, <0.600, 0.741, 0.608>, <0.820, 0.643, 0.537>]
      [0.336, 0.474, <0.820, 0.643, 0.537>, <0.886, 0.780, 0.714>]
      [0.474, 0.810, <0.886, 0.780, 0.714>, <0.996, 0.643, 0.537>]
      [0.810, 0.836, <0.996, 0.643, 0.537>, <0.973, 0.973, 0.976>]
      [0.836, 1.001, <0.973, 0.973, 0.976>, <0.973, 0.973, 0.976>])

//------ Dark Green & Browns
define Grnt27_map
   color_map(
      [0.000, 0.043, <0.773, 0.647, 0.569>, <0.431, 0.322, 0.227>]
      [0.043, 0.113, <0.431, 0.322, 0.227>, <0.278, 0.282, 0.216>]
      [0.113, 0.304, <0.278, 0.282, 0.216>, <0.278, 0.282, 0.216>]
      [0.304, 0.426, <0.278, 0.282, 0.216>, <0.459, 0.341, 0.243>]
      [0.426, 0.843, <0.459, 0.341, 0.243>, <0.459, 0.341, 0.243>]
      [0.843, 0.878, <0.459, 0.341, 0.243>, <0.459, 0.341, 0.243>]
      [0.878, 0.983, <0.459, 0.341, 0.243>, <0.278, 0.282, 0.216>]
      [0.983, 1.001, <0.278, 0.282, 0.216>, <0.773, 0.647, 0.569>])

//------- Red Swirl
define Grnt28_map
   color_map(
      [0.000, 0.155, <0.686, 0.235, 0.282>, <0.686, 0.235, 0.282>]
      [0.155, 0.328, <0.686, 0.235, 0.282>, <0.494, 0.243, 0.294>]
      [0.328, 0.474, <0.494, 0.243, 0.294>, <0.769, 0.329, 0.373>]
      [0.474, 0.647, <0.769, 0.329, 0.373>, <0.769, 0.329, 0.373>]
      [0.647, 0.810, <0.769, 0.329, 0.373>, <0.686, 0.235, 0.282>]
      [0.810, 0.922, <0.686, 0.235, 0.282>, <0.792, 0.388, 0.427>]
      [0.922, 1.001, <0.792, 0.388, 0.427>, <0.686, 0.235, 0.282>])

//-------- White Tan & thin Reds
define Grnt29_map
   color_map(
      [0.000, 0.053, <0.784, 0.627, 0.522>, <0.784, 0.627, 0.624>]
      [0.053, 0.263, <0.784, 0.627, 0.624>, <0.824, 0.557, 0.376>]
      [0.263, 0.281, <0.824, 0.557, 0.376>, <0.643, 0.380, 0.376>]
      [0.281, 0.325, <0.643, 0.380, 0.376>, <0.839, 0.722, 0.722>]
      [0.325, 0.711, <0.839, 0.722, 0.722>, <0.784, 0.627, 0.522>]
      [0.711, 0.798, <0.784, 0.627, 0.522>, <0.769, 0.380, 0.376>]
      [0.798, 0.895, <0.769, 0.380, 0.376>, <0.824, 0.557, 0.376>]
      [0.895, 0.982, <0.824, 0.557, 0.376>, <0.784, 0.627, 0.522>]
      [0.982, 1.001, <0.784, 0.627, 0.522>, <0.784, 0.627, 0.522>])

define Grnt0a_map
      color_map(
	 [0.000, 0.153, <0.729, 0.502, 0.451>, 0.306,
			<0.769, 0.686, 0.592>, 0.792]
	 [0.153, 0.398, <0.769, 0.686, 0.592>,  0.792,
			<0.843, 0.753, 0.718>,  0.396]
	 [0.398, 0.559, <0.843, 0.753, 0.718>,  0.396,
			<0.780, 0.667, 0.561>,  0.976]
	 [0.559, 0.729, <0.780, 0.667, 0.561>,  0.976,
			<0.741, 0.659, 0.576>,  0.820]
	 [0.729, 1.001, <0.741, 0.659, 0.576>, 0.820,
			<0.729, 0.502, 0.451>,  0.306])

//----- Translucent Grnt1
define Grnt1a_map
   color_map(
      [0.000, 0.212, <0.898, 0.898, 0.851>, 0.306, 
		     <0.969, 0.980, 0.875>, 0.792]
      [0.212, 0.424, <0.969, 0.980, 0.875>, 0.792, 
		     <0.859, 0.859, 0.859>, 0.396]
      [0.424, 0.627, <0.859, 0.859, 0.859>, 0.396, 
		     <0.992, 0.922, 0.659>, 0.976]
      [0.627, 0.881, <0.992, 0.922, 0.659>, 0.976, 
		     <0.937, 0.965, 0.902>, 0.820]
      [0.881, 1.001, <0.937, 0.965, 0.902>, 0.820, 
		     <0.898, 0.898, 0.851>, 0.306])

//-----Translucent Grnt2
define Grnt2a_map
   color_map(
      [0.000, 0.144, <0.863, 0.757, 0.596>, 0.596,
		     <0.925, 0.792, 0.714>, 0.349]
      [0.144, 0.288, <0.925, 0.792, 0.714>, 0.349,
		     <0.871, 0.702, 0.659>, 0.784]
      [0.288, 0.644, <0.871, 0.702, 0.659>, 0.784,
		     <0.992, 0.922, 0.659>, 0.498]
      [0.644, 0.983, <0.992, 0.922, 0.659>, 0.498,
		     <0.902, 0.812, 0.714>, 0.722]
      [0.983, 1.001, <0.902, 0.812, 0.714>, 0.722,
		     <0.863, 0.757, 0.596>, 0.596])

//-----Translucent Grnt3
define Grnt3a_map
   color_map(
      [0.000, 0.153, <0.831, 0.631, 0.569>, 0.447,
		     <0.925, 0.831, 0.714>, 0.678]
      [0.153, 0.297, <0.925, 0.831, 0.714>, 0.678,
		     <0.871, 0.702, 0.659>, 0.475]
      [0.297, 0.441, <0.871, 0.702, 0.659>, 0.475,
		     <0.831, 0.631, 0.569>, 0.918]
      [0.441, 0.763, <0.831, 0.631, 0.569>, 0.918,
		     <0.937, 0.882, 0.820>, 0.655]
      [0.763, 1.001, <0.937, 0.882, 0.820>, 0.655,
		     <0.831, 0.631, 0.569>, 0.447])

//-----Translucent Grnt4
define Grnt4a_map
   color_map(
      [0.000, 0.144, <0.804, 0.569, 0.494>, 0.569,
		     <0.816, 0.725, 0.537>, 0.467]
      [0.144, 0.449, <0.816, 0.725, 0.537>, 0.467,
		     <0.820, 0.580, 0.522>, 0.584]
      [0.449, 0.568, <0.820, 0.580, 0.522>, 0.584,
		     <0.882, 0.725, 0.537>, 0.871]
      [0.568, 0.754, <0.882, 0.725, 0.537>, 0.871,
		     <0.855, 0.729, 0.584>, 0.816]
      [0.754, 1.001, <0.855, 0.729, 0.584>, 0.816,
		     <0.804, 0.569, 0.494>, 0.569])

//-----Translucent Grnt4
define Grnt5a_map
   color_map(
      [0.000, 0.178, <0.804, 0.569, 0.494>, 0.569, 
		     <0.855, 0.729, 0.584>, 0.467]
      [0.178, 0.356, <0.855, 0.729, 0.584>, 0.467, 
		     <0.667, 0.502, 0.478>, 0.584]
      [0.356, 0.525, <0.667, 0.502, 0.478>, 0.584, 
		     <0.859, 0.624, 0.545>, 0.871]
      [0.525, 0.729, <0.859, 0.624, 0.545>, 0.871, 
		     <0.855, 0.729, 0.584>, 0.816]
      [0.729, 1.001, <0.855, 0.729, 0.584>, 0.816, 
		     <0.804, 0.569, 0.494>, 0.569])
//-----Translucent Grnt6
define Grnt6a_map
   color_map(
      [0.000, 0.263, <0.545, 0.380, 0.345>, 0.733,
		     <0.588, 0.475, 0.333>, 0.741]
      [0.263, 0.432, <0.588, 0.475, 0.333>, 0.741,
		     <0.675, 0.478, 0.404>, 0.545]
      [0.432, 0.551, <0.675, 0.478, 0.404>, 0.545,
		     <0.757, 0.635, 0.522>, 0.384]
      [0.551, 0.720, <0.757, 0.635, 0.522>, 0.384,
		     <0.659, 0.549, 0.443>, 0.675]
      [0.720, 1.001, <0.659, 0.549, 0.443>, 0.675,
		     <0.545, 0.380, 0.345>, 0.733])
//-----Translucent Grnt7
define Grnt7a_map
   color_map(
      [0.000, 0.119, <0.439, 0.310, 0.282>, 0.631,
		     <0.463, 0.369, 0.259>, 0.847]
      [0.119, 0.322, <0.463, 0.369, 0.259>, 0.847,
		     <0.541, 0.369, 0.298>, 0.549]
      [0.322, 0.449, <0.541, 0.369, 0.298>, 0.549,
		     <0.573, 0.424, 0.286>, 0.965]
      [0.449, 0.729, <0.573, 0.424, 0.286>, 0.965,
		     <0.494, 0.396, 0.306>, 0.741]
      [0.729, 1.001, <0.494, 0.396, 0.306>, 0.741,
		     <0.439, 0.310, 0.282>, 0.631])
//-----Aqua Tints 
define Grnt8a_map
   color_map(
      [0.000, 0.119, <0.310, 0.384, 0.420>, 0.631,
		     <0.322, 0.369, 0.416>, 0.847]
      [0.119, 0.322, <0.322, 0.369, 0.416>, 0.847,
		     <0.424, 0.369, 0.420>, 0.549]
      [0.322, 0.449, <0.424, 0.369, 0.420>, 0.549,
		     <0.373, 0.424, 0.518>, 0.965]
      [0.449, 0.729, <0.373, 0.424, 0.518>, 0.965,
		     <0.482, 0.573, 0.533>, 0.741]
      [0.729, 1.001, <0.482, 0.573, 0.533>, 0.741,
		     <0.310, 0.384, 0.420>, 0.631])
//-----Alpha Creams With Cracks
define Grnt9a_map
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
//--------- Alpha Cream Rose & light yellow
define Grnt10a_map
   color_map(
      [0.000, 0.179, <0.843, 0.655, 0.655>, 0.455,
		     <0.886, 0.769, 0.627>, 0.608]
      [0.179, 0.368, <0.886, 0.769, 0.627>, 0.608,
		     <0.906, 0.820, 0.714>, 0.392]
      [0.368, 0.538, <0.906, 0.820, 0.714>, 0.392,
		     <0.851, 0.671, 0.671>, 0.659]
      [0.538, 0.744, <0.851, 0.671, 0.671>, 0.659,
		     <0.890, 0.792, 0.675>, 0.392]
      [0.744, 0.983, <0.890, 0.792, 0.675>, 0.392,
		     <0.827, 0.612, 0.612>, 0.706]
      [0.983, 1.001, <0.827, 0.612, 0.612>, 0.706,
		     <0.843, 0.655, 0.655>, 0.455])
//--------- Alpha Light Grays
define Grnt11a_map
   color_map(
      [0.000, 0.154, <0.894, 0.886, 0.886>, 0.659,
		     <0.745, 0.745, 0.753>, 0.584]
      [0.154, 0.308, <0.745, 0.745, 0.753>, 0.584,
		     <0.902, 0.902, 0.859>, 0.780]
      [0.308, 0.444, <0.902, 0.902, 0.859>, 0.780,
		     <0.729, 0.706, 0.694>, 0.686]
      [0.444, 0.615, <0.729, 0.706, 0.694>, 0.686,
		     <0.588, 0.592, 0.635>, 0.424]
      [0.615, 0.803, <0.588, 0.592, 0.635>, 0.424,
		     <0.608, 0.616, 0.659>, 0.761]
      [0.803, 1.001, <0.608, 0.616, 0.659>, 0.761,
		     <0.894, 0.886, 0.886>, 0.659])
//--------- Alpha Creams & Tans
define Grnt12a_map
      color_map(
	 [0.000, 0.154, <0.890, 0.690, 0.690>, 0.659,
			<0.996, 0.835, 0.737>, 0.659]
	 [0.154, 0.308, <0.996, 0.835, 0.737>, 0.659,
			<0.745, 0.635, 0.651>, 0.780]
	 [0.308, 0.444, <0.745, 0.635, 0.651>, 0.780,
			<0.733, 0.596, 0.557>, 0.686]
	 [0.444, 0.615, <0.733, 0.596, 0.557>, 0.686,
			<0.996, 0.835, 0.737>, 0.659]
	 [0.615, 0.803, <0.996, 0.835, 0.737>, 0.659,
			<0.765, 0.616, 0.659>, 0.761]
	 [0.803, 1.001, <0.765, 0.616, 0.659>, 0.761,
			<0.890, 0.690, 0.690>, 0.659])
//--------- Alpha Creams & Grays
define Grnt13a_map
   color_map(
      [0.000, 0.154, <0.800, 0.651, 0.557>, 0.000,
		     <0.996, 0.835, 0.737>, 0.608]
      [0.154, 0.308, <0.996, 0.835, 0.737>, 0.608,
		     <0.800, 0.651, 0.557>, 0.635]
      [0.308, 0.444, <0.800, 0.651, 0.557>, 0.635,
		     <0.694, 0.624, 0.604>, 0.294]
      [0.444, 0.615, <0.694, 0.624, 0.604>, 0.294,
		     <0.800, 0.651, 0.557>, 0.816]
      [0.615, 0.812, <0.800, 0.651, 0.557>, 0.816,
		     <0.725, 0.655, 0.651>, 0.957]
      [0.812, 1.001, <0.725, 0.655, 0.651>, 0.957,
		     <0.800, 0.651, 0.557>, 0.000])
//--------- Cream Rose & light moss
define Grnt14a_map
   color_map(
      [0.000, 0.154, <0.996, 0.969, 0.800>, 0.373,
		     <0.996, 0.682, 0.604>, 0.412]
      [0.154, 0.308, <0.996, 0.682, 0.604>, 0.412,
		     <0.906, 0.820, 0.714>, 0.616]
      [0.308, 0.444, <0.906, 0.820, 0.714>, 0.616,
		     <0.816, 0.631, 0.537>, 0.443]
      [0.444, 0.615, <0.816, 0.631, 0.537>, 0.443,
		     <0.890, 0.792, 0.675>, 0.745]
      [0.615, 0.812, <0.890, 0.792, 0.675>, 0.745,
		     <0.973, 0.627, 0.627>, 0.600]
      [0.812, 1.001, <0.973, 0.627, 0.627>, 0.600,
		     <0.996, 0.969, 0.800>, 0.373])
//--------- Alpha Sand & light Orange
define Grnt15a_map
   color_map(
      [0.000, 0.154, <0.996, 0.824, 0.780>, 0.412,
		     <0.996, 0.698, 0.624>, 0.412]
      [0.154, 0.308, <0.996, 0.698, 0.624>, 0.412,
		     <0.906, 0.675, 0.553>, 0.616]
      [0.308, 0.444, <0.906, 0.675, 0.553>, 0.616,
		     <0.996, 0.682, 0.604>, 0.412]
      [0.444, 0.615, <0.996, 0.682, 0.604>, 0.412,
		     <0.996, 0.824, 0.780>, 0.412]
      [0.615, 0.812, <0.996, 0.824, 0.780>, 0.412,
		     <0.973, 0.627, 0.627>, 0.600]
      [0.812, 1.001, <0.973, 0.627, 0.627>, 0.600,
		     <0.996, 0.824, 0.780>, 0.412])
//--------- Cream Rose & light moss
define Grnt16a_map
   color_map(
      [0.000, 0.078, <0.769, 0.722, 0.690>, 0.180,
		     <0.745, 0.690, 0.655>, 1.000]
      [0.078, 0.96, <0.745, 0.690, 0.655>, 1.000,
		     <0.839, 0.804, 0.780>, 1.000]
      [0.96, 1.001, <0.839, 0.804, 0.780>, 0.278,
		     <0.769, 0.722, 0.690>, 0.180])
define Grnt17a_map
   color_map(
      [0.000, 0.034, <0.027, 0.012, 0.012>, 0.000,
		     <0.851, 0.812, 0.741>, 0.235]
      [0.034, 0.342, <0.851, 0.812, 0.741>, 0.235,
		     <0.792, 0.694, 0.690>, 0.839]
      [0.342, 0.462, <0.792, 0.694, 0.690>, 0.839,
		     <0.631, 0.506, 0.471>, 0.608]
      [0.462, 0.632, <0.631, 0.506, 0.471>, 0.608,
		     <0.851, 0.812, 0.741>, 0.922]
      [0.632, 0.983, <0.851, 0.812, 0.741>, 0.922,
		     <0.647, 0.655, 0.655>, 0.282]
      [0.983, 1.001, <0.647, 0.655, 0.655>, 0.282,
		     <0.027, 0.012, 0.012>, 0.000])
define Grnt18a_map
   color_map(
      [0.000, 0.128, <0.820, 0.580, 0.580>, 0.000,
		     <0.851, 0.812, 0.741>, 0.235]
      [0.128, 0.282, <0.851, 0.812, 0.741>, 0.235,
		     <0.792, 0.694, 0.690>, 0.282]
      [0.282, 0.393, <0.792, 0.694, 0.690>, 0.282,
		     <0.647, 0.655, 0.655>, 0.133]
      [0.393, 0.590, <0.647, 0.655, 0.655>, 0.133,
		     <0.851, 0.812, 0.741>, 0.333]
      [0.590, 0.983, <0.851, 0.812, 0.741>, 0.333,
		     <0.647, 0.655, 0.655>, 0.282]
      [0.983, 1.001, <0.647, 0.655, 0.655>, 0.282,
		     <0.820, 0.580, 0.580>, 0.000])
//--------- Gray Marble with White feather Viens with Alpha
define Grnt19a_map
   color_map(
      [0.00, 0.35, White,   0.0, DimGray, 0.5]
      [0.35, 0.40, DimGray, 0.5, DimGray, 0.8]
      [0.40, 1.00, DimGray, 0.8, DimGray, 0.9])
//--------- White Feature Viens 
define Grnt20a_map
   color_map(
      [0.0, 0.3, white, 0.0, white, 0.7]
      [0.3, 0.4, white, 0.7, white, 1.0]
      [0.4, 1.0, white, 1.0, white, 1.0])
//--------- Thinner White Feature Viens 
define Grnt21a_map
   color_map(
      [0.0, 0.2, white, 0.4, white, 0.8]
      [0.2, 0.3, white, 0.8, white, 1.0]
      [0.3, 1.0, white, 1.0, white, 1.0])
define Grnt22a_map
   color_map(
      [0.000, 0.175, <1.000, 0.718, 0.541>, 0.890,
		     <0.843, 0.678, 0.655>, 0.753]
      [0.175, 0.228, <0.843, 0.678, 0.655>, 0.753,
		     <0.906, 0.831, 0.773>, 0.98]
      [0.228, 0.386, <0.906, 0.831, 0.773>, 0.698,
		     <0.992, 0.718, 0.545>, 0.794]
      [0.386, 0.412, <0.992, 0.718, 0.545>, 0.794,
		     <0.333, 0.188, 0.067>, 0.784]
      [0.412, 0.439, <0.333, 0.188, 0.067>, 0.784,
		     <0.925, 0.557, 0.514>, 0.778]
      [0.439, 0.684, <0.925, 0.557, 0.514>, 0.678,
		     <0.984, 0.780, 0.655>, 0.696]
      [0.684, 0.781, <0.984, 0.780, 0.655>, 0.696,
		     <0.965, 0.847, 0.675>, 0.880]
      [0.781, 0.982, <0.965, 0.847, 0.675>, 0.880,
		     <1.000, 0.718, 0.541>, 0.990]
      [0.982, 1.001, <1.000, 0.718, 0.541>, 0.890,
		     <1.000, 0.718, 0.541>, 0.890])
//---------- Transparent Green Moss Colors
define Grnt23a_map
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
define Grnt24a_map
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

//--------------------Crack & OverTint /Red---------------------
define Crack1_map
   color_map(
      [0.00, 0.03, white,     0.6, white,   0.6]
      [0.03, 0.30, <0.6, 0.1, 0.1>, 0.8, <0.6, 0.1, 0.1>, 0.8] 
      [0.30, 0.95, <0.6, 0.1, 0.1>, 0.8, dimgray, 0.8] 
      [0.95, 1.00, dimgray,   0.8, white,   1.0])

//-----Alpha  Dark Cracks
define Crack2_map
   color_map(
      [0.00, 0.50, white, 1.0, white, 1.0]
      [0.50, 0.54, white, 1.0, white, 0.0]
      [0.54, 1.00, white, 1.0, white, 1.0])

//---------- Overtint Green with Black Cracks
define Crack3_map
   color_map(
      [0.00, 0.10, white,     0.60, white,     1.0]
      [0.10, 0.90, DarkGreen, 0.80, DarkGreen, 0.90] 
      [0.90, 1.00, white,     0.90, white,     1.0])

//--------- Overtint with White Crack
define Crack4_map
   color_map(
      [0.00, 0.03, white, 0.30, white,   1.0]
      [0.03, 0.97, white, 0.70, DimGray, 0.9] 
      [0.97, 1.00, white, 0.90, white,   1.0])

define Grnt0 Grnt0_map[granite_fn_04]
define Grnt1 Grnt1_map[granite_fn_06]
define Grnt2 Grnt2_map[granite_fn_05]
define Grnt3 Grnt3_map[granite_fn_05]
define Grnt4 Grnt4_map[granite_fn_05]
define Grnt5 Grnt5_map[granite_fn_05]
define Grnt6 Grnt6_map[granite_fn_05]
define Grnt7 Grnt7_map[granite_fn_05]
define Grnt8 Grnt8_map[granite_fn_06]
define Grnt9 Grnt9_map[granite_fn_06]
define Grnt10 Grnt10_map[granite_fn_06]
define Grnt11 Grnt11_map[granite_fn_06]
define Grnt12 Grnt12_map[granite_fn_06]
define Grnt13 Grnt13_map[granite_fn_06]
define Grnt14 Grnt14_map[granite_fn_06]
define Grnt15 Grnt15_map[granite_fn_06]
define Grnt16 Grnt16_map[granite_fn_06]
define Grnt17 Grnt17_map[granite_fn_06]
define Grnt18 Grnt18_map[granite_fn_04]
define Grnt19 Grnt19_map[granite_fn]
define Grnt20 Grnt20_map[granite_fn]
define Grnt21 Grnt21_map[granite_fn]
define Grnt22 Grnt22_map[granite_fn]
define Grnt23 Grnt23_map[granite_fn_05]
define Grnt24 Grnt24_map[granite_fn_06]
define Grnt25 Grnt25_map[marble_fn]
define Grnt26 Grnt26_map[granite_fn]
define Grnt27 Grnt27_map[granite_fn]
define Grnt28 Grnt28_map[granite_fn_07]
define Grnt29 Grnt29_map[noisy_layers]

define Grnt0a Grnt0a_map[granite_fn_06]
define Grnt1a Grnt1a_map[granite_fn_06]
define Grnt2a Grnt2a_map[granite_fn_06]
define Grnt3a Grnt3a_map[granite_fn_06]
define Grnt4a Grnt4a_map[granite_fn_06]
define Grnt5a Grnt5a_map[granite_fn_05]
define Grnt6a Grnt6a_map[granite_fn_06]
define Grnt7a Grnt7a_map[granite_fn_06]
define Grnt8a Grnt8a_map[granite_fn_06]
define Grnt9a Grnt9a_map[granite_fn_06]
define Grnt10a Grnt10a_map[granite_fn_06]
define Grnt11a Grnt11a_map[granite_fn_06]
define Grnt12a Grnt12a_map[granite_fn_06]
define Grnt13a Grnt13a_map[granite_fn_06]
define Grnt14a Grnt14a_map[granite_fn_06]
define Grnt15a Grnt15a_map[granite_fn_06]
define Grnt16a Grnt16a_map[granite_fn_06]
define Grnt17a Grnt17a_map[granite_fn_06]
define Grnt18a Grnt18a_map[granite_fn_06]
define Grnt19a Grnt19a_map[granite_fn_06]
define Grnt20a Grnt20a_map[granite_fn]
define Grnt21a Grnt21a_map[granite_fn]
define Grnt22a Grnt22a_map[granite_fn_05]
define Grnt23a Grnt23a_map[granite_fn]
define Grnt24a Grnt24a_map[granite_fn_05]

define Crack1
   texture {
      special shiny {
	 color Crack1_map[noisy_layers2]
	 diffuse 0.4
	 // Filter lower layers by crack color
	 transmission Crack1_map[noisy_layers2], 1, 1
	 }
      }

define Crack2
   texture {
      special shiny {
	 color Crack2_map[noisy_layers2]
	 transmission Crack2_map[noisy_layers2], 1, 1
	 }
      }
define Crack3 Crack3_map[noisy_layers]
define Crack4 Crack4_map[noisy_layers]

//------------- Deep Rose & Green Marble with large White Swirls
define Stone1
texture {
   layered
      texture { special shiny { color Grnt0a }
	        scale <2, 3, 2> rotate <0, 0, -30> },
      texture { special shiny { color Grnt7 }
	        scale <2, 3, 2> rotate <0, 0, 40> }
   }

//------------- Light Greenish Tan Marble with Agate style veining
define Stone2
texture {
   layered
      texture { special shiny { color Grnt7a }
	        scale <1.3, 2, 1.3> rotate <0, 0, -30> rotate <40, 0, 0> },
      texture { special shiny { color Grnt0 }
	        scale <2, 3, 2> rotate <0, 0, 40> }
      }

//------------- Rose & Yellow Marble with fog white veining
define Stone3
texture {
   layered
      texture { special shiny { color Grnt0a }
	        scale <2, 3, 2> rotate <0, 0, -30> },
      texture { special shiny { color Grnt5 }
	        scale <2, 3, 2> rotate <0, 0, 40> }
   }

//------------- Tan Marble with Rose patches
define Stone4
texture {
   layered
      texture { special shiny { color Grnt10a }
		scale <1, 3, 2> rotate <0, 0, -30> },
      texture { special shiny { color Grnt6 }
		scale <1.5, 3, 2> rotate <0, 0, 40> }
   }

//------------- White Cream Marble with Pink veining
define Stone5
texture {
   layered
      Crack1 { scale <1, 2, 1.4> rotate <10, 0, -20> },
      texture { special shiny { color Grnt17a }
	        scale <2, 3, 2> rotate <0, 0, -30> },
      texture { special shiny { color Grnt12 }
		scale <2, 3, 2> rotate <0, 0, 40> }
   }

//------------- Rose & Yellow Cream Marble
define Stone6
texture {
   layered
      Crack1 { scale <1, 2, 1.4> rotate <10, 0, -20> },
      texture { special shiny { color Grnt19a }
	        scale <2, 4, 1.3> rotate <0, 0, 30> },
      texture { special shiny { color Grnt18 }
	        scale <1.5, 3, 3> rotate <0, 0, 40> }
   }

//------------- Light Coffee Marble with darker patches
define Stone7
texture {
   layered
      texture { special shiny { color Grnt9a }
	        scale <3.5, 5, 4> rotate <0, 0, 60> },
      texture { special shiny { color Grnt6a }
	        scale <1, 3, 2> rotate <0, 0, 40> },
      texture { shiny { color salmon } }
   }

//------------- Gray Granite with white patches
define Stone8
texture {
   layered
      texture { special shiny { color Grnt9a }
	        scale <5, 3, 4> rotate <0, 0, 40> },
      texture { special shiny { color Grnt0a }
	        scale <2, 3, 2> rotate <0, 0, -30> },
      shiny_white
   }

//------------- White & Light Blue Marble with light violets
define Stone9
texture {
   layered
      Crack1 { scale <1, 2, 1.4> rotate <10, 0, -20> },
      texture { special shiny { color Grnt9 }
	        scale <1.2, 2.4, 1.2> rotate <0, 0, -30> rotate <40, 0, 0> }
   }

//------------- Dark Brown & Tan swirl Granite with gray undertones
define Stone10
texture {
   layered
      texture { special shiny { color Grnt3a }
	        scale <1, 2, 1> rotate <0, 0, -50> },
      texture { special shiny { color Grnt17a }
	        scale <3, 6, 2> rotate <0, 0, 50> },
      shiny_black
   }

//------------- Rose & White Marble with dark tan swirl
define Stone11
texture {
   layered
      Crack1 { scale <1, 2, 1.4> rotate <10, 0, -20> },
      texture { special shiny { color Grnt2a }
	        scale <3, 3, 4> rotate <0, 0, 40> },
      texture { special shiny { color Grnt15a }
	        scale <1.2, 3, 1.5> rotate <70, 0, 30> },
      shiny_black
   }

//------------- White & Pinkish Tan Marble
define Stone12
texture {
   layered
      texture { special shiny { color Grnt0a }
	        scale <1, 3, 2> rotate <0, 0, -30> },
      texture { special shiny { color Grnt23 }
	        scale <1, 5, 1> rotate <0, 0, 50> }
   }

//------------- Medium Gray Blue Marble
define Stone13
texture {
   layered
      texture { special shiny { color Grnt8a }
	        scale <1, 3, 2> rotate <0, 0, -30> },
      texture { special shiny { color Grnt24 }
	        scale <2, 5, 2> rotate <0, 0, 50> }
   }

//------------- Tan & Olive Marble with gray white veins
define Stone14
texture {
   layered
      texture { special shiny { color Grnt19a }
	        scale <1, 3, 1> rotate <0, 0, 40> },
      texture { special shiny { color Grnt6 }
	        scale <2, 3, 2> rotate <0, 0, -30> }
   }

//------------- Deep Gray Marble with white veining
define Stone15
texture {
   layered
      texture { special shiny { color Grnt8a }
	        scale <1, 2, 1> rotate <0, 0, 40> },
      texture { special shiny { color Grnt20 }
	        scale <1, 2, 2> rotate <0, 0, -30> }
   }

//------------- Peach & Yellow Marble with white veining
define Stone16
texture {
   layered
      texture { special shiny { color Grnt20a }
	        scale <1, 2, 1> rotate <0, 0, 40> },
      texture { special shiny { color Grnt19a }
	        scale <2, 4, 2> rotate <0, 0, -30> },
      texture { special shiny { color Grnt18 }
	        scale <1.3, 2, 2> rotate <0, 0, -30> }
   }

//------------- White Marble with gray veining
define Stone17
texture {
   layered
      texture { special shiny { color Crack3 }
	        scale <1, 2, 1.4> rotate <10, 0, -20> },
      texture { special shiny { color Grnt20 }
	        scale <1, 2, 2> rotate <0, 0, -30> }
   }

//------------- Green Jade with white veining
define Stone18
texture {
   layered
      texture { special shiny { color Crack4 }
	        scale <0.7, 1, 1> rotate <10, 0, -20> },
      texture { special shiny { color Grnt20a }
	        scale <2.5, 2, 0.5> rotate <0, 0, -50> },
      texture { special shiny { color Grnt22 }
	        scale <1.5, 0.7, 0.5> rotate <0, 0, 40> },
      texture { special shiny { color SeaGreen } }
   }

//------------- Peach Granite with white patches & green trim
define Stone19
texture {
   layered
      texture { special shiny { color Grnt20a }
	        scale <2, 3, 1> rotate <10, 0, -20> },
      texture { special shiny { color Grnt26 }
	        scale <1, 0.7, 0.5> rotate <0, 0, 40> }
   }

//------------- Brown & Olive Marble with white veining
define Stone20
texture {
   layered
      Crack1 { scale <1, 0.6, 1> rotate <10, 0, -20> },
      texture { special shiny { color Grnt20a }
	        scale <1.9, 3, 0.5> rotate <0, 0, -50> },
      texture { special shiny { color Grnt12a }
	        scale <1, 1.3, 2> rotate <0, 0, 40> },
      texture { special shiny { color Grnt27 }
	        scale <0.7, 0.99, 0.7> rotate <0, 0, 40> }
   }

//------------- Red Marble with gray & white veining
define Stone21
texture {
   layered
      texture { special shiny { color Crack4 }
	        scale <1, 0.6, 1> rotate <10, 0, -20> },
      texture { special shiny { color Grnt22 }
	        scale <1, 2, 2> rotate <0, 0, 40> },
      texture { special shiny { color Grnt28 }
	        scale <1.3, 2.5, 1.7> rotate <0, 0, 40> }
   }

//------------- Dark Tan Marble with gray & white veining
define Stone22
texture {
   layered
      texture { special shiny { color Crack4 }
	        rotate <10, 0, -40> },
      texture { special shiny { color Grnt22 }
	        scale <2, 4, 1.5> rotate <0, 0, -50> },
      texture { special shiny { color Grnt8a }
	        scale <1, 2, 2> rotate <0, 0, 40> },
      texture { shiny { color tan } }
   }

//------------- Peach & Cream Marble with orange veining
define Stone23
texture {
   layered
      clear_yellow,
      Crack1 { scale <1, 2, 1.5> rotate <0, 0, 40> },
      texture { special shiny { color Grnt24a }
	        scale <2, 1, 2> rotate <40, 0, 30> },
      texture { special shiny { color Grnt29 }
	        scale <1, 1, 2> rotate <40, 0, 30> }
   }

//------------- Green & Tan Moss Marble
define Stone24
texture {
   layered
      texture { special shiny { color Grnt23a }
	        scale <2, 1, 2> rotate <40, 0, -30> },
      texture { special shiny { color Grnt25 }
	        scale <1, 1, 2> rotate <0, 0, 50> rotate <20, 0, 30> }
   }

