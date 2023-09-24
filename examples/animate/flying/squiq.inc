define rcs_spike_1_10 object { sphere <0, 0, 1>, 1 scale <1, 1, 10> }
define rcs_spike_3_20 object { sphere <0, 0, 1>, 1 scale <3, 3, 20> }
define rcs_spike_5_15 object { sphere <0, 0, 1>, 1 scale <5, 5, 15> }

define pie_slice
   object {
      disc <0, 0, 0>, <0, 1, 0>, 20, 40
      uv_bounds 0, 1/4, 0, 1
      rotate <0,-45, 0>
      texture {
	 noise surface {
	    position_fn 3
	    ambient 0.4
	    diffuse 0.4
	    specular white, 0.5
	    microfacet Reitz 5
	    color_map(
	       [20, 25, black, 1, blue, 0]
	       [25, 30, blue, yellow]
	       [30, 35, yellow, magenta]
	       [35, 40, magenta, 0, black, 1])
	    }
	 rotate <90, 0, 0>
	 }
      }

define bumpy_green
texture {
   special surface {
      normal N + 2 * (dnoise(3*W) - <0.5, 0.5, 0.5>)
      ambient green, 0.1
      diffuse green, 0.6
      specular white, 0.5
      microfacet Cook 5
      }
   scale <0.1, 0.1, 0.1>
   }

define plane_blue
texture {
   surface {
      ambient <0.14,0.42,0.56>, 0.2
      diffuse <0.14,0.42,0.56>, 0.6
      specular white, 0.5
      reflection white, 0.3
      microfacet Reitz 20
      }
   }

define black_yellow_stripes
texture {
   noise surface {
      ambient 0.1
      diffuse 0.6
      specular white, 0.5
      microfacet Cook 5
      normal 1
      position_fn 1
      lookup_fn 1
      color_map(
	 [0, 0.5, black, black]
	 [0.5, 1, yellow, yellow])
      }
   scale <0.1, 0.1, 0.1>
   rotate <0, 90, 0>
   }

define template1
object {
    object { sphere <0, 0, 0>, 5 plane_blue }
   + rcs_spike_3_20 { rotate <0, 30, 0> }
   + rcs_spike_3_20 { rotate <0,-30, 0> }
   + rcs_spike_1_10 { translate <0, 0, 20> rotate 0.5 * <0, 30, 0> }
   + rcs_spike_1_10 { translate <0, 0, 20> rotate 0.5 * <0,-30, 0> }
   + rcs_spike_1_10 { translate <0, 0, 20> rotate 1.5 * <0, 30, 0> }
   + rcs_spike_1_10 { translate <0, 0, 20> rotate 1.5 * <0,-30, 0> }
   + rcs_spike_5_15 { rotate <0, 180, 0> black_yellow_stripes }
   + pie_slice
   bumpy_green
   }
