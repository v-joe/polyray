//
// Build the segment counter "segment" and linear interpolating variable
// "increment" for each of the component pieces of the animation.
//

// Define how many frames are used in each segment of the animation
define segment_lengths [200, 50, 50, 100]

// Define arrays to hold the start/stop frames of each segment
define start_times
   [0, segment_lengths[0], segment_lengths[0] + segment_lengths[1],
    segment_lengths[0] + segment_lengths[1] + segment_lengths[2]]
define end_times   [start_times[0]+segment_lengths[0],
		    start_times[1]+segment_lengths[1],
		    start_times[2]+segment_lengths[2],
		    start_times[3]+segment_lengths[3]]

// Tell Polyray how many frames will be rendered
start_frame  start_times[0]
total_frames segment_lengths[0] + segment_lengths[1] + segment_lengths[2] +
	     segment_lengths[3]
end_frame total_frames

// Determine which sequence we are doing
define segment (frame < start_times[1] ? 0 :
		  (frame < start_times[2] ? 1 :
		     (frame < start_times[3] ? 2 : 3)))


// Define the linear interpolation variable "increment" for this segment
define t0 start_times[segment]
define t1 end_times[segment]
define increment (frame - t0) / (t1 - t0)
