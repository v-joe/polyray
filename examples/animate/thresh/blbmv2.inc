if (frame == start_frame) {
   // Set the initial position and velocity of each blob component
   static define px0  0
   static define py0  0
   static define vx0 0.5 * brownian(<0, 0, 0>, 1)[0]
   static define vy0 0.5 * brownian(<0, 0, 0>, 1)[1]

   static define px1  0
   static define py1  0
   static define vx1 0.5 * brownian(<0, 0, 0>, 1)[0]
   static define vy1 0.5 * brownian(<0, 0, 0>, 1)[1]

   static define px2  0
   static define py2  0
   static define vx2 0.5 * brownian(<0, 0, 0>, 1)[0]
   static define vy2 0.5 * brownian(<0, 0, 0>, 1)[1]

   static define px3  0
   static define py3  0
   static define vx3 0.5 * brownian(<0, 0, 0>, 1)[0]
   static define vy3 0.5 * brownian(<0, 0, 0>, 1)[1]

   static define px4  0
   static define py4  0
   static define vx4 0.5 * brownian(<0, 0, 0>, 1)[0]
   static define vy4 0.5 * brownian(<0, 0, 0>, 1)[1]

   static define px5  0
   static define py5  0
   static define vx5 0.5 * brownian(<0, 0, 0>, 1)[0]
   static define vy5 0.5 * brownian(<0, 0, 0>, 1)[1]

   static define px6  0
   static define py6  0
   static define vx6 0.5 * brownian(<0, 0, 0>, 1)[0]
   static define vy6 0.5 * brownian(<0, 0, 0>, 1)[1]

   static define px7  0
   static define py7  0
   static define vx7 0.5 * brownian(<0, 0, 0>, 1)[0]
   static define vy7 0.5 * brownian(<0, 0, 0>, 1)[1]

   static define px8  0
   static define py8  0
   static define vx8 0.5 * brownian(<0, 0, 0>, 1)[0]
   static define vy8 0.5 * brownian(<0, 0, 0>, 1)[1]

   static define px9  0
   static define py9  0
   static define vx9 0.5 * brownian(<0, 0, 0>, 1)[0]
   static define vy9 0.5 * brownian(<0, 0, 0>, 1)[1]
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

   static define px5 px5 + vx5
   if (px5 > bwidth) {
      static define px5 2 * bwidth - px5
      static define vx5 -vx5 + rbounce
      }
   else if (px5 < -bwidth)  {
      static define px5 -2 * bwidth - px5
      static define vx5 -vx5 + rbounce
      }
   static define py5 py5 + vy5
   if (py5 > bheight) {
      static define py5 2 * bheight - py5
      static define vy5 max(min(-vy5 + rbounce, 1), -1)
      }
   else if (py5 < -bheight) {
      static define py5 -2 * bheight - py5
      static define vy5 max(min(-vy5 + rbounce, 1), -1)
      }

   static define px6 px6 + vx6
   if (px6 > bwidth) {
      static define px6 2 * bwidth - px6
      static define vx6 -vx6 + rbounce
      }
   else if (px6 < -bwidth)  {
      static define px6 -2 * bwidth - px6
      static define vx6 -vx6 + rbounce
      }
   static define py6 py6 + vy6
   if (py6 > bheight) {
      static define py6 2 * bheight - py6
      static define vy6 max(min(-vy6 + rbounce, 1), -1)
      }
   else if (py6 < -bheight) {
      static define py6 -2 * bheight - py6
      static define vy6 max(min(-vy6 + rbounce, 1), -1)
      }

   static define px7 px7 + vx7
   if (px7 > bwidth) {
      static define px7 2 * bwidth - px7
      static define vx7 -vx7 + rbounce
      }
   else if (px7 < -bwidth)  {
      static define px7 -2 * bwidth - px7
      static define vx7 -vx7 + rbounce
      }
   static define py7 py7 + vy7
   if (py7 > bheight) {
      static define py7 2 * bheight - py7
      static define vy7 max(min(-vy7 + rbounce, 1), -1)
      }
   else if (py7 < -bheight) {
      static define py7 -2 * bheight - py7
      static define vy7 max(min(-vy7 + rbounce, 1), -1)
      }

   static define px8 px8 + vx8
   if (px8 > bwidth) {
      static define px8 2 * bwidth - px8
      static define vx8 -vx8 + rbounce
      }
   else if (px8 < -bwidth)  {
      static define px8 -2 * bwidth - px8
      static define vx8 -vx8 + rbounce
      }
   static define py8 py8 + vy8
   if (py8 > bheight) {
      static define py8 2 * bheight - py8
      static define vy8 max(min(-vy8 + rbounce, 1), -1)
      }
   else if (py8 < -bheight) {
      static define py8 -2 * bheight - py8
      static define vy8 max(min(-vy8 + rbounce, 1), -1)
      }

   static define px9 px9 + vx9
   if (px9 > bwidth) {
      static define px9 2 * bwidth - px9
      static define vx9 -vx9 + rbounce
      }
   else if (px9 < -bwidth)  {
      static define px9 -2 * bwidth - px9
      static define vx9 -vx9 + rbounce
      }
   static define py9 py9 + vy9
   if (py9 > bheight) {
      static define py9 2 * bheight - py9
      static define vy9 max(min(-vy9 + rbounce, 1), -1)
      }
   else if (py9 < -bheight) {
      static define py9 -2 * bheight - py9
      static define vy9 max(min(-vy9 + rbounce, 1), -1)
      }
   }
