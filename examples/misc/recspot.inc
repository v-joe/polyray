// Lighting for the recursive room

define from0 < 4, 0.5, 7.3>  define at0 < 2, 4.5, 10>
define from1 <-4, 0.5, 7.3>  define at1 <-2, 4.5, 10>
define from2 < 5.3, 0.5, 1>  define at2 <8, 4.5, 3>
define from3 < 5.3, 0.5, 9>  define at3 <8, 4.5, 6>
define from4 <-5.3, 0.5, 1>  define at4 <-8, 4.5, 3>
define from5 <-5.3, 0.5, 9>  define at5 <-8, 4.5, 6>

define from6 < 7.5, 9.5, 9>   define at6 <-1, 0,-1>
define from7 <-7.5, 9.5, 9>   define at7 < 1, 0,-1>
define from8 < 7.5, 9.5, 0>   define at8 <-1, 0, 1>
define from9 <-7.5, 9.5, 0>   define at9 < 1, 0, 1>

define dif0 (at0 - from0)/|at0 - from0|
define dif1 (at1 - from1)/|at1 - from1|
define dif2 (at2 - from2)/|at2 - from2|
define dif3 (at3 - from3)/|at3 - from3|
define dif4 (at4 - from4)/|at4 - from4|
define dif5 (at5 - from5)/|at5 - from5|
define dif6 (at6 - from6)/|at6 - from6|
define dif7 (at7 - from7)/|at7 - from7|
define dif8 (at8 - from8)/|at8 - from8|
define dif9 (at9 - from9)/|at9 - from9|

define spotlight_texture
texture {
   special surface {
      normal N + dnoise(3*W)
      ambient grey, 0.2
      diffuse grey, 0.3
      specular white, 0.6
      microfacet Cook 10
      }
   scale <0.05, 0.05, 0.05>
   }

define bright_white texture { surface { color white ambient 1 diffuse 0 } }

// Spots pointed towards the back wall
define s0 0.5
define s1 0.3
define s2 0.6
define l0 1.0
define l1 0.3
spot_light white, from0, at0, 1, 15, 25
spot_light white, from1, at1, 1, 15, 25
object {
   (object { sphere from0, s0 } & object { disc from0+l1*dif0, dif0, s0 })
   + object { cone from0, s1, from0 + l0 * dif0, s2 }
   + object { disc from0+l0*dif0, dif0, s0 bright_white shading_flags 8 }
   spotlight_texture
   }
object {
   (object { sphere from1, s0 } & object { disc from1+l1*dif1, dif1, s0 })
   + object { cone from1, s1, from1 + l0 * dif1, s2 }
   + object { disc from1+l0*dif1, dif1, s0 bright_white shading_flags 8 }
   spotlight_texture
   }

// Spots pointed towards the right wall
spot_light white, from2, at2, 1, 15, 25
spot_light white, from3, at3, 1, 15, 25
object {
   (object { sphere from2, s0 } & object { disc from2+l1*dif2, dif2, s0 })
   + object { cone from2, s1, from2 + l0 * dif2, s2 }
   + object { disc from2+l0*dif2, dif2, s0 bright_white shading_flags 8 }
   spotlight_texture
   }
object {
   (object { sphere from3, s0 } & object { disc from3+l1*dif3, dif3, s0 })
   + object { cone from3, s1, from3 + l0 * dif3, s2 }
   + object { disc from3+l0*dif3, dif3, s0 bright_white shading_flags 8 }
   spotlight_texture
   }

// Spots pointed towards the left wall
spot_light white, from4, at4, 1, 15, 25
spot_light white, from5, at5, 1, 15, 25
object {
   (object { sphere from4, s0 } & object { disc from4+l1*dif4, dif4, s0 })
   + object { cone from4, s1, from4 + l0 * dif4, s2 }
   + object { disc from4+l0*dif4, dif4, s0 bright_white shading_flags 8 }
   spotlight_texture
   }
object {
   (object { sphere from5, s0 } & object { disc from5+l1*dif5, dif5, s0 })
   + object { cone from5, s1, from5 + l0 * dif5, s2 }
   + object { disc from5+l0*dif5, dif5, s0 bright_white shading_flags 8 }
   spotlight_texture
   }

spot_light white, from6, at6, 1, 5, 7
spot_light white, from7, at7, 1, 5, 7
spot_light white, from8, at8, 1, 5, 7
spot_light white, from9, at9, 1, 5, 7
object {
   (object { sphere from6, s0 } & object { disc from6+l1*dif6, dif6, s0 })
   + object { cone from6, s1, from6 + l0 * dif6, s0 }
   + object { disc from6+l0*dif6, dif6, s0 bright_white shading_flags 8 }
   spotlight_texture
   }
object {
   (object { sphere from7, s0 } & object { disc from7+l1*dif7, dif7, s0 })
   + object { cone from7, s1, from7 + l0 * dif7, s0 }
   + object { disc from7+l0*dif7, dif7, s0 bright_white shading_flags 8 }
   spotlight_texture
   }
