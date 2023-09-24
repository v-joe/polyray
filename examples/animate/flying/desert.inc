//
// Add a floor
//
define ground_ripple
texture {
   noise surface {
      color <1, 0.66, 0.2>
      normal 2
      frequency 100
      bump_scale 2
      ambient 0.3
      diffuse 0.7
      }
   }
object {
   disc <0, 0, 0>, <0, 1, 0>, 10001
   ground_ripple { scale <200, 200, 200> }
   }

//
// Things to fly around
//
object {
   cylinder <-8, 0, 0>, <-8, 7, 0>, 0.5
   texture { checker matte_white, matte_black }
   }
object { sphere <-8, 7, 0>, 1 mirror }
object {
   cylinder <8, 0, 0>, <8, 7, 0>, 0.5
   texture { checker matte_white, matte_black }
   }
object { sphere <8, 7, 0>, 1 mirror }
