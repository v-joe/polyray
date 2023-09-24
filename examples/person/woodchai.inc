// Let's make a chair, extends from <0, 0, 0> to <3, 6, 3>
define leg_wood reddish_wood { rotate <90, 0, 0> }
define wooden_chair
object {
   // Legs
   object { box <0, 0,   0>, <0.4, 6,   0.4> leg_wood } +
   object { box <0, 0,   0>, <0.4, 2.5, 0.4> translate <2.2, 0, 0  > leg_wood} +
   object { box <0, 0,   0>, <0.4, 2.5, 0.4> translate <2.2, 0, 2.8> leg_wood} +
   object { box <0, 0,   0>, <0.4, 6,   0.4> translate <0,   0, 2.8> leg_wood} +
   // Top bars
   object { box <0, 5.6, 0>, <0.4, 6,   3.2> reddish_wood } +
   object { box <0, 4.8, 0>, <0.4, 5.2, 3.2> reddish_wood } +
   object { box <0, 4,   0>, <0.4, 4.4, 3.2> reddish_wood } +
   object { box <0, 3.2, 0>, <0.4, 3.6, 3.2> reddish_wood } +
   // Seat
   object {
      box <0, 2.5, 0>, <2.6, 2.7, 3.2>
      reddish_wood {
	 rotate <90, 0, 0>
	 translate <1.3, 0, 1.6>
	 }
      }
   }
