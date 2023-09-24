define wafer_height 0.4

define wafer0_width peg_base_width/1.1
define wafer1_width wafer0_width * 0.8
define wafer2_width wafer1_width * 0.8

define wafer0
object {
    object { disc <0, 0,   0>, <0, -1, 0>, peg_width, wafer0_width }
  + object { disc <0, wafer_height, 0>, <0,  1, 0>, peg_width, wafer0_width }
  + object { cylinder <0, 0, 0>, <0, wafer_height, 0>, wafer0_width }
  + object { cylinder <0, 0, 0>, <0, wafer_height, 0>, peg_width }
  jade
  }

define wafer1
object {
    object { disc <0, 0,   0>, <0, -1, 0>, peg_width, wafer1_width }
  + object { disc <0, wafer_height, 0>, <0,  1, 0>, peg_width, wafer1_width }
  + object { cylinder <0, 0, 0>, <0, wafer_height, 0>, wafer1_width }
  + object { cylinder <0, 0, 0>, <0, wafer_height, 0>, peg_width }
  jade { translate <30, 30, 30> }
  }

define wafer2
object {
    object { disc <0, 0,   0>, <0, -1, 0>, peg_width, wafer2_width }
  + object { disc <0, wafer_height, 0>, <0,  1, 0>, peg_width, wafer2_width }
  + object { cylinder <0, 0, 0>, <0, wafer_height, 0>, wafer2_width }
  + object { cylinder <0, 0, 0>, <0, wafer_height, 0>, peg_width }
  jade { translate <20, 17, 20> }
  }

// Which wafer are we working on this move?
define working_wafer [2, 1, 2, 0, 2, 1, 2, 2]

// Define how high the wafer is when starting the move
// and when ending the move
define wafer_positions [[[-1, 0], [-1, 1], [-1, 2]],
                        [[-1, 0], [-1, 1], [ 0, 0]],
                        [[-1, 0], [ 1, 0], [ 0, 0]],
                        [[-1, 0], [ 1, 0], [ 1, 1]],
                        [[ 0, 0], [ 1, 0], [ 1, 1]],
                        [[ 0, 0], [ 1, 0], [-1, 0]],
                        [[ 0, 0], [ 0, 1], [-1, 0]],
                        [[ 0, 0], [ 0, 1], [ 0, 2]]]

if (submove == 0) {
   // We are moving the wafer from the start position, straight up
   define moving_wafer_height
      peg_base_height +
      (1 - submove_increment) * wafer_height * 
         wafer_positions[move][working_wafer[move]][1] +
      submove_increment * peg_length * 1.1
   define moving_wafer_offset
      peg_offset * wafer_positions[move][working_wafer[move]][0]
   }
else if (submove == 1) {
   // We are moving the wafer from the start peg over to the next peg
   define moving_wafer_height peg_base_height + peg_length * 1.1
   define moving_wafer_offset
      peg_offset *
         ((1 - submove_increment) *
                    wafer_positions[move][working_wafer[move]][0] +
          submove_increment * wafer_positions[move+1][working_wafer[move]][0])
   }
else if (submove == 2) {
   // We are moving the wafer from the top of the end peg straight down
   define moving_wafer_height
      peg_base_height +
      submove_increment * wafer_height *
         wafer_positions[move+1][working_wafer[move]][1] +
      (1 - submove_increment) * peg_length * 1.1
   define moving_wafer_offset
      peg_offset * wafer_positions[move+1][working_wafer[move]][0]
   }
else if (submove == 3) {
   // We are moving only the arm from the end position, straight up
   define moving_wafer_height
      peg_base_height +
      (1 - submove_increment) * wafer_height * 
         wafer_positions[move+1][working_wafer[move]][1] +
      submove_increment * peg_length * 1.1
   define moving_wafer_offset
      peg_offset * wafer_positions[move+1][working_wafer[move]][0]
   }
else if (submove == 4) {
   // We are moving only the arm from the end peg over to the new start peg
   define moving_wafer_height peg_base_height + peg_length * 1.1
   define moving_wafer_offset
      peg_offset *
         ((1 - submove_increment) *
                    wafer_positions[move+1][working_wafer[move]][0] +
          submove_increment * wafer_positions[move][working_wafer[move+1]][0])
   }
else if (submove == 5) {
   // We are moving only the arm from the top of the new start peg straight down
   define moving_wafer_height
      peg_base_height +
      submove_increment * wafer_height *
         wafer_positions[move][working_wafer[move+1]][1] +
      (1 - submove_increment) * peg_length * 1.1
   define moving_wafer_offset
      peg_offset * wafer_positions[move][working_wafer[move+1]][0]
   }

define which_wafer_width0 ((working_wafer[move] == 0) ? wafer0_width :
			   ((working_wafer[move] == 1) ? wafer1_width :
			    wafer2_width))
define which_wafer_width1 ((working_wafer[move+1] == 0) ? wafer0_width :
			   ((working_wafer[move+1] == 1) ? wafer1_width :
			    wafer2_width))

if (submove < 3)
   define wafer_width which_wafer_width0
else
   define wafer_width which_wafer_width1

if (working_wafer[move] == 0) {
   if (submove < 3) {
      define wafer_height0 moving_wafer_height
      define wafer_offset0 moving_wafer_offset
      }
   else {
      define wafer_height0 peg_base_height +
			   wafer_height * wafer_positions[move+1][0][1]
      define wafer_offset0 peg_offset * wafer_positions[move+1][0][0]
      }
   define wafer_height1 peg_base_height +
                        wafer_height * wafer_positions[move][1][1]
   define wafer_offset1 peg_offset * wafer_positions[move][1][0]
   define wafer_height2 peg_base_height +
                        wafer_height * wafer_positions[move][2][1]
   define wafer_offset2 peg_offset * wafer_positions[move][2][0]
   }
else if (working_wafer[move] == 1) {
   define wafer_height0 peg_base_height +
                        wafer_height * wafer_positions[move][0][1]
   define wafer_offset0 peg_offset * wafer_positions[move][0][0]
   if (submove < 3) {
      define wafer_height1 moving_wafer_height
      define wafer_offset1 moving_wafer_offset
      }
   else {
      define wafer_height1 peg_base_height +
			   wafer_height * wafer_positions[move+1][1][1]
      define wafer_offset1 peg_offset * wafer_positions[move+1][1][0]
      }
   define wafer_height2 peg_base_height +
                        wafer_height * wafer_positions[move][2][1]
   define wafer_offset2 peg_offset * wafer_positions[move][2][0]
   }
else {
   define wafer_height0 peg_base_height +
                        wafer_height * wafer_positions[move][0][1]
   define wafer_offset0 peg_offset * wafer_positions[move][0][0]
   define wafer_height1 peg_base_height +
                        wafer_height * wafer_positions[move][1][1]
   define wafer_offset1 peg_offset * wafer_positions[move][1][0]
   if (submove < 3) {
      define wafer_height2 moving_wafer_height
      define wafer_offset2 moving_wafer_offset
      }
   else {
      define wafer_height2 peg_base_height +
			   wafer_height * wafer_positions[move+1][2][1]
      define wafer_offset2 peg_offset * wafer_positions[move+1][2][0]
      }
   }
