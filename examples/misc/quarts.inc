// Bicorn
define a 1.0
define bicorn
object {
    object { polynomial  y^2*(a^2 - (x^2 + z^2)) - (x^2 + z^2 + 2*a*y - a^2)^2 }
  & object { box <-3, -3, -0.3>, <3, 3, 0.3> }
  bounding_box <-3, -3, -0.3>, <3, 3, 0.3>
  }

// Bifolia
define a 3
define bifolia
object {
     object { polynomial (x^2 + y^2 + z^2)^2 - a*(x^2 + z^2)*y }
   & object { box <-1, -1, -0.1>, <1, 1, 1> }
   bounding_box <-1, -1, -0.1>, <1, 1, 1>
   }

// Ovals of Cassini rotated around the y axis
define a 0.45
define b 0.5
define c 16
define cassini
object {
     object { polynomial (x^2 + y^2 + z^2 + a^2)^2 - c*a^2*(x^2 + z^2) - b^2 }
   & object { box <-2, -2, -0.6>, <2, 2, 2> }
   bounding_box <-2, -2, -0.6>, <2, 2, 2>
   scale <2/3, 2/3, 2/3>
   }

// A cubic saddle.
define csaddle
object {
     object { polynomial x^3 - y^3 - z }
   & object { box <-1, -1, -1>, <1, 1, 1> }
   bounding_box <-1, -1, -1>, <1, 1, 1>
   scale <5/6, 5/6, 5/6>
   }

// Devils curve variant
define devils_curve
object {
     object { polynomial x^4 + 2*x^2*z^2 - 0.36*x^2 - y^4 + 0.25*y^2 + z^4 }
   & object { box <-2, -2, -0.5>, <2, 2, 0.5> }
   bounding_box <-2, -2, -0.5>, <2, 2, 0.5>
   rotate <90, 0, 0>
   scale <2/3, 2/3, 2/3>
   }

// Folium
define a 1
define b 1
define folium
object {
     object { polynomial (y^2 + z^2) * (1 + (b - 4*a)*x) + x^2*(1 + b)
              scale <0.1, 0.1, 0.1> }
   & object { box <-1, -1, -1>, <1, 1, 1> }
   bounding_box <-1, -1, -1>, <1, 1, 1>
   }

// Hyperbolic torus
define r0 0.6
define r1 0.4
define hyp_torus
object {
   object { polynomial x^4 + 2*x^2*y^2 - 2*x^2*z^2 - 2*(r0^2+r1^2)*x^2 + y^4 -
              2*y^2*z^2 + 2*(r0^2-r1^2)*y^2 + z^4 + 2*(r0^2+r1^2)*z^2 +
              (r0^2-r1^2)^2 }
   &  object { box <-2, -2, -2>, <2, 2, 2> }
   bounding_box <-2, -2, -2>, <2, 2, 2>
   scale <2/3, 2/3, 2/3>
   }

// Kampyle of Eudoxus
define a 0.2
define c 1
define kampyle
object {
   object { polynomial (y^2 + z^2) - c^2 * x^4 + c^2 * a^2 * x^2 }
   & object { box <-1, -1, -1>, <1, 1, 1> }
   bounding_box <-1, -1, -1>, <1, 1, 1>
   }

// Lemniscate of Gerono
define lemniscate
object {
   polynomial x^4 - x^2 + y^2 + z^2
   bounding_box <-1, -1, -1>, <1, 1, 1>
   }

// Two parallel sheets with a blob on the top
define a 10
define c 1
define qloop
object {
   object { polynomial (x^2+y^2+z^2+a*c*x)^2-(x^2+y^2+z^2)*(c-a*x)^2 }
   & object { box <-2, -1, -0.4>, <2, 1, 0.4> }
   bounding_box <-2, -1, -0.4>, <2, 1, 0.4>
   rotate <0, 0, -90>
   rotate <90, 0, 0>
   }

// Monkey Saddle
define c 1
define monkey_saddle
object {
   object { polynomial c*(x^3 - 3*x*y^2) - z }
   & object { box <-1, -1, -1>, <1, 1, 1> }
   bounding_box <-1, -1, -1>, <1, 1, 1>
   }

// Parabolic Torus
define r0 0.6
define r1 0.5
define par_torus
object {
   object {
   polynomial x^4 + 2*x^2*y^2 - 2*x^2*z - (r0^2+r1^2)*x^2 + y^4 - 2*y^2*z +
              (r0^2-r1^2)*y^2 + z^2 + (r0^2+r1^2)*z + (r0^2-r1^2)^2 }
   & object { box <-2, -2, -2>, <2, 2, 2> }
   bounding_box <-2, -2, -2>, <2, 2, 2>
   scale <2/3, 2/3, 2/3>
   }

// Piriform - very nice teardrop shape
define piriform
object {
   object { polynomial (x^4 - x^3) + y^2 + z^2 translate <-0.5, 0, 0> }
   bounding_box <-0.6, -0.6, -0.6>, <0.6, 0.6, 0.6>
   scale <2, 1, 1>
   rotate <0, 90, 0>
   translate <0, 1, 0>
   }

// Quartic paraboloid.  Looks like the quadric paraboloid, but is squared
// off on the bottom and sides.
define quart_parab
object {
   object { polynomial x^4 + z^4 - y }
   & object { box <-2, 0, -0.6>, <2, 2, 2> }
   bounding_box <-2, 0, -0.6>, <2, 2, 2>
   scale <2/3, 2/3, 2/3>
   rotate <90, 0, 0>
   }

// Quartic saddle - looks like the quadric saddle, but is squared off in
// the middle.  The equation is:
define quart_saddle
object {
   object { polynomial x^4 - z^4 - y }
   & object { box <-1.5, -1.5, -1>, <1.5, 1.3, 1.5> }
   bounding_box <-1.5, -1.5, -1>, <1.5, 1.3, 1.5>
   scale <1/2, 1/2, 1/2>
   rotate <90, 0, 0>
   }

// Space Needle?
define a 1
define c 0.1
define space_needle
object {
   object { polynomial (x^2 + z^2) * y^2 + c^2 * (x^2 + z^2) - c^2 * a^2 }
   & object { box <-1, -3, -1>, <1, 3, 1> }
   bounding_box <-1, -3, -1>, <1, 3, 1>
   scale <1/2, 1/2, 1/2>
   rotate <90, 0, 0>
   }

// Quartic parabola - a 4th degree polynomial (has two bumps at the bottom)
// that has been swept around the z axis.
define parab4
object {
   object { polynomial 0.1*(x^2 + z^2)^2 - (x^2 + z^2) + 0.9 - y
            scale <0.5, 0.5, 0.5> }
   & object { box <-2, -2, -0.7>, <2, 0.7, 2> }
   bounding_box <-2, -2, -0.7>, <2, 0.7, 2>
   scale <2/3, 2/3, 2/3>
   rotate <90, 0, 0>
   }

// Steiners quartic surface.
define steiner
object {
   object { polynomial x^2*y^2 + x^2*z^2 + y^2*z^2 + x*y*z }
   bounding_box <-sqrt(2), -sqrt(2), -sqrt(2)>, <sqrt(2), sqrt(2), sqrt(2)>
   scale <2, 2, 2>
   rotate <-40, -10, 0>
   }

// Strophoid
define a 1
define b -0.1
define c 0.4
define strophoid
object {
   object { polynomial (b - x)*(y^2 + z^2) - c^2*a*x^2 - c^2*x^3 }
   & object { box <-1, -1, -1>, <1, 1, 1> }
   bounding_box <-1, -1, -1>, <1, 1, 1>
   rotate <0, 0, -90>
   }

// A drop coming out of water? This is a curve formed by using the equation
// x = 1/2 y^2 (y + 1) as the radius of a cylinder having the y-axis as
// its central axis.
define tear3
object {
   object { polynomial y^2 + z^2 - 0.5*(x^3 + x^2) }
   & object { box <-1, -1, -1>, <1, 1, 1> }
   bounding_box <-1, -1, -1>, <1, 1, 1>
   }

// Glob - sort of like basic teardrop shape.
define tear5
object {
   object { polynomial  0.5*x^5 + 0.5*x^4 - (y^2 + z^2) }
   & object { box <-1, -1, -1>, <1, 1, 1> }
   bounding_box <-1, -1, -1>, <1, 1, 1>
   }

// Torus - basic doughnut shape.  The distance from the origin to the center
// of the ring is "r0", the distance from the center of the ring to the surface
// is "r1".
define r0 0.5
define r1 0.2
define torus_5_2
object {
   object { polynomial (x^2 + y^2 + z^2 - (r0^2 + r1^2))^2 - 4*r0^2*(r1^2-z^2) }
   bounding_box <-(r0+r1), -(r0+r1), -r1>, <r0+r1, r0+r1, r1>
   rotate <90, 0, 0>
   }

// Crossed Trough This is a surface with four pieces that sweep up from
// the x-z plane.
define crossed_trough
object {
   object { polynomial x^2*z^2 - y }
   & object { box <-2, 0, -2>, <2, 2, 2> }
   bounding_box <-2, 0, -2>, <2, 2, 2>
   scale <1/2, 1/2, 1/2>
   }

// Variant of a lemniscate - the two lobes are much more teardrop-like.
define a 1
define c 2
define twin_glob
object {
   object { polynomial y^2 + z^2 - c^2*a^2*x^4 + c^2*x^6 }
   & object { box <-1, -1, -1>, <1, 1, 1> }
   bounding_box <-1, -1, -1>, <1, 1, 1>
   }

// Witch of Agnesi
define a 0.04
define witch_agnesi
object {
   object { polynomial a * (y - 1) + (x^2 + z^2) * y }
   & object { sphere <0, 0, 0>, 1 }
   bounding_box <-1, -1, -1>, <1, 1, 1>
   }
