define shiny_red
texture {
   surface {
      ambient red, 0.2
      diffuse red, 0.8
      specular white, 0.7
      microfacet Reitz 10
      }
   }

object { disc <0, -1, 0>, <0, 1, 0>, 30 }
object {
   sphere <0, 0, 0>, 1
   shiny_red
   uv_steps 32, 64
   }
