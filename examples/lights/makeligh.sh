polyray alight1.pi   -o alight1.tga
mv out000.tga alig0.tga
mv out001.tga alig1.tga
mv out002.tga alig2.tga
polyray alight2.pi   -o alight2.tga
polyray flare1.pi    -o flare1.tga
polyray flare4.pi
mv out040.tga flare4.tga
polyray depspot.pi -r 1
mv out001.tga depspot.tga
rm out000.tga
