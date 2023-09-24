// Define the set of surfaces that range from opaque blue to transparent
// red.  These are used as a particle ages to make it seem to fade away
if (frame == start_frame) {
   static define surf0
   surface {
      ambient 1.0 * 0.1
      diffuse 1.0 * 0.6
      specular white, 0.0 * 0.6
      transmission white, 0.0, 1
      }

   static define surf1
   surface {
      ambient 0.9 * 0.1
      diffuse 0.9 * 0.6
      specular white, 0.9 * 0.6
      transmission white, 0.1, 1
      }

   static define surf2
   surface {
      ambient 0.8 * 0.1
      diffuse 0.8 * 0.6
      specular white, 0.8 * 0.6
      transmission white, 0.2, 1
      }

   static define surf3
   surface {
      ambient 0.7 * 0.1
      diffuse 0.7 * 0.6
      specular white, 0.7 * 0.6
      transmission white, 0.3, 1
      }

   static define surf4
   surface {
      ambient 0.6 * 0.1
      diffuse 0.6 * 0.6
      specular white, 0.6 * 0.6
      transmission white, 0.4, 1
      }

   static define surf5
   surface {
      ambient 0.5 * 0.1
      diffuse 0.5 * 0.6
      specular white, 0.5 * 0.6
      transmission white, 0.5, 1
      }

   static define surf6
   surface {
      ambient 0.4 * 0.1
      diffuse 0.4 * 0.6
      specular white, 0.4 * 0.6
      transmission white, 0.6, 1
      }

   static define surf7
   surface {
      ambient 0.3 * 0.1
      diffuse 0.3 * 0.6
      specular white, 0.3 * 0.6
      transmission white, 0.7, 1
      }

   static define surf8
   surface {
      ambient 0.2 * 0.1
      diffuse 0.2 * 0.6
      specular white, 0.2 * 0.6
      transmission white, 0.8, 1
      }

   static define surf9
   surface {
      ambient 0.1 * 0.1
      diffuse 0.1 * 0.6
      specular white, 0.1 * 0.6
      transmission white, 0.9, 1
      }

   static define tex0 texture { surf0 { color <0.00, 0, 1.00> } }
   static define tex1 texture { surf1 { color <0.10, 0, 0.90> } }
   static define tex2 texture { surf2 { color <0.20, 0, 0.80> } }
   static define tex3 texture { surf3 { color <0.30, 0, 0.70> } }
   static define tex4 texture { surf4 { color <0.40, 0, 0.60> } }
   static define tex5 texture { surf5 { color <0.50, 0, 0.50> } }
   static define tex6 texture { surf6 { color <0.60, 0, 0.40> } }
   static define tex7 texture { surf7 { color <0.70, 0, 0.30> } }
   static define tex8 texture { surf8 { color <0.80, 0, 0.20> } }
   static define tex9 texture { surf9 { color <0.90, 0, 0.10> } }

   static define bsphere0 object { sphere <0, 0, 0>, 0.20 tex0 uv_steps 6, 3 }
   static define bsphere1 object { sphere <0, 0, 0>, 0.25 tex1 uv_steps 6, 3 }
   static define bsphere2 object { sphere <0, 0, 0>, 0.30 tex2 uv_steps 6, 3 }
   static define bsphere3 object { sphere <0, 0, 0>, 0.35 tex3 uv_steps 6, 3 }
   static define bsphere4 object { sphere <0, 0, 0>, 0.40 tex4 uv_steps 6, 3 }
   static define bsphere5 object { sphere <0, 0, 0>, 0.45 tex5 uv_steps 6, 3 }
   static define bsphere6 object { sphere <0, 0, 0>, 0.50 tex6 uv_steps 6, 3 }
   static define bsphere7 object { sphere <0, 0, 0>, 0.55 tex7 uv_steps 6, 3 }
   static define bsphere8 object { sphere <0, 0, 0>, 0.60 tex8 uv_steps 6, 3 }
   static define bsphere9 object { sphere <0, 0, 0>, 0.65 tex9 uv_steps 6, 3 }
   }
