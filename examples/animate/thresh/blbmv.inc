if (frame == start_frame) {
   // Set the initial position and velocity of each blob component
   static define px0  0
   static define py0  0
   static define vx0  0.25
   static define vy0 -0.15

   static define px1  0
   static define py1  0
   static define vx1 -0.05
   static define vy1  0.3

   static define px2  0
   static define py2  0
   static define vx2  0.25
   static define vy2  0.15

   static define px3  0
   static define py3  0
   static define vx3  0.2
   static define vy3 -0.15

   static define px4  0
   static define py4  0
   static define vx4 -0.1
   static define vy4  0.2
   }
else {
   // Recalculate the positions and velocities of each component

   static define px0 px0 + vx0
   if (px0 > bwidth) {
      static define px0 2 * bwidth - px0
      static define vx0 -vx0 + rbounce
      }
   else if (px0 < -bwidth)  {
      static define px0 -2 * bwidth - px0
      static define vx0 -vx0 + rbounce
      }
   static define py0 py0 + vy0
   if (py0 > bheight) {
      static define py0 2 * bheight - py0
      static define vy0 max(min(-vy0 + rbounce, 1), -1)
      }
   else if (py0 < -bheight) {
      static define py0 -2 * bheight - py0
      static define vy0 max(min(-vy0 + rbounce, 1), -1)
      }

   static define px1 px1 + vx1
   if (px1 > bwidth) {
      static define px1 2 * bwidth - px1
      static define vx1 -vx1 + rbounce
      }
   else if (px1 < -bwidth)  {
      static define px1 -2 * bwidth - px1
      static define vx1 -vx1 + rbounce
      }
   static define py1 py1 + vy1
   if (py1 > bheight) {
      static define py1 2 * bheight - py1
      static define vy1 max(min(-vy1 + rbounce, 1), -1)
      }
   else if (py1 < -bheight) {
      static define py1 -2 * bheight - py1
      static define vy1 max(min(-vy1 + rbounce, 1), -1)
      }

   static define px2 px2 + vx2
   if (px2 > bwidth) {
      static define px2 2 * bwidth - px2
      static define vx2 -vx2 + rbounce
      }
   else if (px2 < -bwidth)  {
      static define px2 -2 * bwidth - px2
      static define vx2 -vx2 + rbounce
      }
   static define py2 py2 + vy2
   if (py2 > bheight) {
      static define py2 2 * bheight - py2
      static define vy2 max(min(-vy2 + rbounce, 1), -1)
      }
   else if (py2 < -bheight) {
      static define py2 -2 * bheight - py2
      static define vy2 max(min(-vy2 + rbounce, 1), -1)
      }

   static define px3 px3 + vx3
   if (px3 > bwidth) {
      static define px3 2 * bwidth - px3
      static define vx3 -vx3 + rbounce
      }
   else if (px3 < -bwidth)  {
      static define px3 -2 * bwidth - px3
      static define vx3 -vx3 + rbounce
      }
   static define py3 py3 + vy3
   if (py3 > bheight) {
      static define py3 2 * bheight - py3
      static define vy3 max(min(-vy3 + rbounce, 1), -1)
      }
   else if (py3 < -bheight) {
      static define py3 -2 * bheight - py3
      static define vy3 max(min(-vy3 + rbounce, 1), -1)
      }

   static define px4 px4 + vx4
   if (px4 > bwidth) {
      static define px4 2 * bwidth - px4
      static define vx4 -vx4 + rbounce
      }
   else if (px4 < -bwidth)  {
      static define px4 -2 * bwidth - px4
      static define vx4 -vx4 + rbounce
      }
   static define py4 py4 + vy4
   if (py4 > bheight) {
      static define py4 2 * bheight - py4
      static define vy4 max(min(-vy4 + rbounce, 1), -1)
      }
   else if (py4 < -bheight) {
      static define py4 -2 * bheight - py4
      static define vy4 max(min(-vy4 + rbounce, 1), -1)
      }
   }
