// The location of the chamber of the cannon
static define CPos <-10, 0, -8>

// Locations of the three cannon balls - they all start inside the cannon
static define Ball1_Pos CPos
static define Ball2_Pos CPos
static define Ball3_Pos CPos
static define Ball4_Pos CPos
static define Ball5_Pos CPos
static define Ball6_Pos CPos

static define ang1 1 * pi / 12 // 15 degrees
static define ang2 2 * pi / 12 // 30 
static define ang3 3 * pi / 12 // 45
static define ang4 4 * pi / 12 // 60
static define ang5 5 * pi / 12 // 75

static define Ball1_Vel <cos(ang4)*cos(ang4),sin(ang4),sin(ang4)*cos(ang4)>
static define Ball2_Vel <cos(ang2)*cos(ang4),sin(ang4),sin(ang2)*cos(ang4)>
static define Ball3_Vel <cos(ang3)*cos(ang3),sin(ang3),sin(ang3)*cos(ang3)>
static define Ball4_Vel <cos(ang5)*cos(ang2),sin(ang2),sin(ang5)*cos(ang2)>
static define Ball5_Vel <cos(ang1)*cos(ang3),sin(ang3),sin(ang1)*cos(ang3)>
static define Ball6_Vel <cos(ang3)*cos(ang2),sin(ang2),sin(ang3)*cos(ang2)>

static define Gravity <0, -0.1, 0>

static define bumpy_gray
texture {
   special surface {
      color grey
      normal N + (dnoise(3*W) - white/2)
      ambient 0.2
      diffuse 0.3
      specular white, 0.7
      microfacet Cook 5
      }
   scale <0.05, 0.05, 0.05>
   }

static define Cannon
object {
   (  object { sphere <0, 0, 0>, 1.75 }
    + object { cylinder <0, 0, 0>, <4, 0, 0>, 1.3 }
    + object { torus 1.3, 0.3, <4, 0, 0>, <1, 0, 0> })
   - object { cylinder <-0.5, 0, 0>, <4.5, 0, 0>, 1 }
   texture { shiny { color <0.14, 0.42, 0.56> } }
   }

static define Projectile
object {
   sphere <0, 0, 0>, 1
   bumpy_gray
   }
