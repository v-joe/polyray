#
# Note that you need a C compiler to use a few of the examples here.
# The Symantec C/C++ compiler is used here, you will need to change
# the command line to whatever you happen to own.
#
polyray hfnoise.pi  -o hfnoise.tga
polyray sinfn.pi    -o sinfn.tga
polyray sombfn.pi   -o sombfn.tga
sc -ml -Id:/sc/include wake.c
wake
polyray wake.pi     -o wake.tga
sc -ml -Id:/sc/include islgen.c
islgen
polyray islands.pi  -o islands.tga -x 256 -y 256
sc -ml -Id:/sc/include spiral.c
spiral
polyray spiral.pi   -o spiral.tga -x 256 -y 256
polyray circ.pi     -o circ.tga -M 512
copy maze2 maze2.tga
polyray maze.pi     -o maze.tga
rm maze2.tga
polyray allgrey.pi  -o allgrey.tga -p 8
polyray marbback.pi -o marbback.tga -p 8
polyray sphland1.pi -o sphland1.tga -p 8
polyray sphland.pi  -o sphland.tga -p 8
polyray shtst1.pi   -o shtst1.tga
polyray shtst2.pi   -o shtst2.tga
polyray shtst3.pi   -o shtst3.tga
polyray shtst4.pi   -o shtst4.tga
polyray shtst5.pi   -o shtst5.tga
convert out*.tga /c5
mv anim.flc shtst5.flc
rm out*.tga
polyray shtst6.pi   -o shtst6.tga
polyray shtst7.pi   -o shtst7.tga
polyray shtst8.pi   -o shtst8.tga
polyray shtst9.pi   -o shtst9.tga
polyray shtst10.pi  -o shtst10.tga
polyray chtst1.pi   -o chtst1.tga
polyray chtst2.pi   -o chtst2.tga
polyray chtst3.pi   -o chtst3.tga
polyray chtst7.pi   -o chtst7.tga
polyray chtst8.pi   -o chtst8.tga
polyray chtst10.pi  -o chtst10.tga
