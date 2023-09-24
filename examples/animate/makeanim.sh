#
# Build the frames of some simple animations.
#
# In order to build the many frames for the animations you will need
# approximately half a Megabyte for each of the animations.  To collect
# the frames into an animation, you will need DTA. You will also need
# a player such as play or dfv to show the resulting animation.
#
cd spline
polyray spline1.pi
convert out*.tga spline1.gif
rm out*.tga
polyray spline2.pi
convert out*.tga spline2.gif
rm out*.tga
cd ../environ
# Must set the quality flag to allow two-sided surfaces in order to see
# the effects of the light on the central sphere of 'ilight.pi'
polyray ilight.pi -q 63
convert out*.tga ilight.gif
rm out*.tga
polyray environ.pi
convert out*.tga environ.gif
rm out*.tga
cd ../particle
polyray part1.pi
convert out*.tga part1.gif
rm out*.tga
polyray part2.pi
convert out*.tga part2.gif
rm out*.tga
polyray part3.pi
convert out*.tga part3.gif
rm out*.tga
polyray part6.pi
convert out*.tga part6.gif
rm out*.tga
cd ../cannon
polyray cannon.pi
convert can*.tga cannon.gif
rm can*.tga
cd ../displace
polyray disp1.pi -r 1 -q 63 -x 200 -y 200
convert out*.tga disp1.gif
rm out*.tga
polyray legen.pi -r 1 -q 40
convert out*.tga legen.gif
rm out*.tga
polyray sphdsp3.pi -r 1
convert out*.tga sphdsp.gif
rm out*.tga
cd ../whirl
polyray whirl.pi
convert whirl*.tga whirl1.gif
rm *.tga
polyray busyb1.pi -r 1 -q 0
convert out*.tga whirl2.gif
rm out*.tga
polyray texpts.pi
convert out*.tga texpts.gif
rm out*.tga
cd ../plane
polyray plane.pi
convert plane*.tga plane.gif
rm *.tga
cd ../texanim
polyray texanim.pi
convert out*.tga texanim.gif
rm out*.tga
polyray wireanim.pi
convert out*.tga wireanim.gif
rm out*.tga
polyray phaserb.pi
convert out*.tga phaserb.gif
rm out*.tga
polyray flame.pi
convert out*.tga flame.gif
rm out*.tga
cd ../ufo
polyray ufo1.pi
convert out*.tga ufo1.gif
rm out*.tga
polyray ufo2.pi
convert out*.tga ufo2.gif
rm out*.tga
polyray ufolaunc.pi
convert out*.tga ufolaunc.gif
rm out*.tga
cd ../street
polyray street.pi
convert out*.tga street.gif
rm out*.tga
cd ../bezier
polyray bezier2.pi -r 1
convert out*.tga bez2.gif
rm out*.tga
cd ../thresh
polyray squirm.pi
convert sqrm*.tga sqrm.gif
rm sqrm*.tga
polyray movesph4.pi
convert out*.tga mvsph4.gif
rm out*.tga
polyray mirblb.pi
convert mirb*.tga mirblb.gif
rm mirb*.tga
polyray posneg.pi -a 2
convert posng*.tga posng.gif
rm posng*.tga
cd ../flying
polyray noiseanm.pi
convert out*.tga noiseanm.gif
rm out*.tga
polyray flyplane.pi
convert out*.tga flyplane.gif
rm out*.tga
cd ../robarm
polyray robarm.pi
convert out*.tga robarm.gif
rm out*.tga
cd ../screw
polyray screw1.pi
convert out*.tga screw.gif
rm out*.tga
