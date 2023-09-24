polyray csg1.pi     -o csg1.tga
polyray csg2.pi     -o csg2.tga
polyray csg3.pi     -o csg3.tga
polyray csg2.pi -r 7 > temp.inc
polyray doraw.pi    -o doraw.tga
rm temp.inc
