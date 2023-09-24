define lego_2_2
object {
     object { box <0, 0, 0>, <2, 1, 2> }
   + object { cylinder <0.5, 1, 0.5>, <0.5, 1.2, 0.5>, 0.4 }
   + object { cylinder <1.5, 1, 0.5>, <1.5, 1.2, 0.5>, 0.4 }
   + object { cylinder <0.5, 1, 1.5>, <0.5, 1.2, 1.5>, 0.4 }
   + object { cylinder <1.5, 1, 1.5>, <1.5, 1.2, 1.5>, 0.4 }
   + object { disc <0.5, 1.2, 0.5>, <0, 1, 0>, 0.4 }
   + object { disc <1.5, 1.2, 0.5>, <0, 1, 0>, 0.4 }
   + object { disc <0.5, 1.2, 1.5>, <0, 1, 0>, 0.4 }
   + object { disc <1.5, 1.2, 1.5>, <0, 1, 0>, 0.4 }
   }

define lego_3_2
object {
     object { box <0, 0, 0>, <3, 1, 2> }
   + object { cylinder <0.5, 1, 0.5>, <0.5, 1.2, 0.5>, 0.4 }
   + object { cylinder <1.5, 1, 0.5>, <1.5, 1.2, 0.5>, 0.4 }
   + object { cylinder <2.5, 1, 0.5>, <2.5, 1.2, 0.5>, 0.4 }
   + object { cylinder <0.5, 1, 1.5>, <0.5, 1.2, 1.5>, 0.4 }
   + object { cylinder <1.5, 1, 1.5>, <1.5, 1.2, 1.5>, 0.4 }
   + object { cylinder <2.5, 1, 1.5>, <2.5, 1.2, 1.5>, 0.4 }
   + object { disc <0.5, 1.2, 0.5>, <0, 1, 0>, 0.4 }
   + object { disc <1.5, 1.2, 0.5>, <0, 1, 0>, 0.4 }
   + object { disc <2.5, 1.2, 0.5>, <0, 1, 0>, 0.4 }
   + object { disc <0.5, 1.2, 1.5>, <0, 1, 0>, 0.4 }
   + object { disc <1.5, 1.2, 1.5>, <0, 1, 0>, 0.4 }
   + object { disc <2.5, 1.2, 1.5>, <0, 1, 0>, 0.4 }
   }

define lego_4_2
object {
   lego_2_2 + lego_2_2 { translate <2, 0, 0> }
   }
