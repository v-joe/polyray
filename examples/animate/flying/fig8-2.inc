define start_times [0, total_frames/4, total_frames/2, 3*total_frames/4]
define end_times [total_frames/4, total_frames/2,
		  3*total_frames/4, total_frames]
define segment floor(4 * frame / total_frames)

define t0 start_times[segment]
define t1 end_times[segment]
//
// Figure out how far along we are in the animation
//
define increment (frame - t0) / (t1 - t0)
//
// Define some values for the movement
//
define start_positions [<-10, 5, 0>, <0, 5, 0>, <10, 5, 0>, <0, 5, 0>]
define end_positions [<0, 5, 0>, <10, 5, 0>, <0, 5, 0>, <-10, 5, 0>]

define start_velocities [<0, 0, 10>, <20, 0, -20>, <0, 0, 10>, <-20, 0, -20>]
define end_velocities [<20, 0, -20>, <0, 0, 10>, <-20, 0, -20>, <0, 0, 10>]

define start_banks [-40, 0, 40, 0]
define end_banks [0, 40, 0, -40]

define pos0 start_positions[segment]
define pos1 end_positions[segment]
define vel0 start_velocities[segment]
define vel1 end_velocities[segment]
define bank0 start_banks[segment]
define bank1 end_banks[segment]
//
// Calculate the various coefficients of the spline path
//
define u  increment
define u2 (u * u)
define u3 (u2 * u)
define a0 (2 * u3 - 3 * u2 + 1)
define a1 (-2 * u3 + 3 * u2)
define a2 (u3 - 2 * u2 + u)
define a3 (u3 - u2)
//
// Finally the position is calculated from the coefficients:
//
define pos (a0 * pos0 + a1 * pos1 + a2 * vel0 + a3 * vel1)
//
// Now move the object
//
define increment_vel (a0 * vel0 + a1 * vel1)
define norm_vel increment_vel / fabs(increment_vel)
define half_angle (norm_vel + <1, 0, 0>) / 2
if (|half_angle| < 0.001)
   define half_angle <0, 0, 1>

define bank_angle (a0 * bank0 + a1 * bank1)
