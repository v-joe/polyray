//
// Determine the camera direction and angle.  What we look at is based on
// the individual segments of the animation.
//
if (segment == 0) {
   // Counterclockwise pan.  Starts by looking up and away from the
   // location of the launch pad.  At the end of segment 0 we are looking
   // down towards the pad, centered on the hole that the UFO will come out of.
   define view_from  <0, 8, -15>
   define view_ang 30
   define view_at view_from + rotate(<0, 0,-15>, <0, increment * -180, 0>) +
		              increment * <0, -6, 0>
   }
else if (segment == 1) {
   // The UFO rises up from the center hole in the launch pad.  We don't do
   // any camera work during this segment.
   define view_from  <0, 8, -15>
   define view_ang 30
   define view_at <0, 2, 0>
   }
else if (segment == 2) {
   // As the UFO finishes it's climb, we keep the camera pointed directly
   // at it.  The field of view is dropped from 30 degrees at the start of
   // this segment down to 15 degrees at the end.
   define view_from  <0, 8, -15>
   define view_ang 30 - increment * 15
   define view_at <0, 2, 0> + increment * <0, 3, 0>
   }
else {
   // This segment simply keeps the camera a constant distance behind the
   // UFO as it flys across the landscape.  The value "100" below is dependent
   // on the values used for the UFO start_pos and end_pos values for this
   // segment.
   define view_from <0, 8, -15+increment*100>
   define view_ang 15
   define view_at <0, 5, increment*100>
   }

// Build the camera according to the values generated above
viewpoint {
   from  view_from
   at view_at
   up <0, 1, 0>
   angle view_ang
   resolution 320, 200
   aspect 4/3
   }
