// Polyray include file for person.pi
// This file contains the rotations for all the joints of the body,
// and color of the figure.
// All rotations follow the left hand rule.
// +x is to the right, +y is to the rear, and +z is up.
// Created by Will Wagner

define headrot		<20,0,0>		// <-45/45,-10/10,-80/80>
define neckrot		<10,0,0>		// <-30/30,-10/10,-20/20>
define torsorot		<0,0,0>		// <-10/100,-25/25,-15/15>
define abdomenrot	<0,0,0>		// <0,-5/5,0>

// Left arm rotations
define laupperrot	<0,30,0>		// <-180/45,-45/180,-90/45>
define lalowerrot	<0,-90,0>		// <-165/0,0,0>
define lhandrot		<0,45,0>		// <-20/20,-90/20,-90/180>
define lfingerrot	<0,0,0>		// <0,-90/0,0>

// Right arm rotations
define raupperrot	<0,-30,0>		// <-180/45,-180/45,-45/90>
define ralowerrot	<0,90,0>		// <-165/0,0,0>
define rhandrot		<0,-45,0>		// <-20/20,-20/90,-180/90>
define rfingerrot	<0,0,0>		// <0,0/90,0>

// Left leg rotations
define lthighrot	<0,0,0>		// <-90/15,0/45,-45/0>
define lcalfrot		<0,0,0>		// <0/160,0,0>
define lfootrot		<0,0,0>		// <-5/50,-5/5,-10/10>

// Right leg rotations
define rthighrot	<0,0,0>		// <-90/15,-45/0,0/45>
define rcalfrot		<0,0,0>		// <0/160,0,0>
define rfootrot		<0,0,0>		// <-5/50,-5/5,-10/10>
