#
# Batch file to process all of the Polyray data files.
# 
# Processing all of the data files can take one or more days, depending
# on the speed of your processor.  If you want a quick look, then
# move to a subdirectory and look for an appropriate test image to render.
#
# You need to have POLYRAY_PATH set to a directory that contains colors.inc
# and texture.inc.  For example you could have this in autoexec.bat:
#
#    set POLYRAY_PATH=d:/polyray/dat
#
# or you could simply uncomment the line above in this batch file (after
# fixing the directory reference).
#
# In a couple of cases, a C compiler is necessary to generate data,
# in particular for the height field.  The batch files are currently
# set to use the Symantec C++ compiler, however a simple change
# of the command line should be all you need for another compiler.
#
# In order to generated the animations you will need the program 'dta' on
# your path.
#
cd simple
sh makesimp.sh
cd ../lights
sh makeligh.sh
cd ../flare
sh makeflar.sh
cd ../raw
sh makeraw.sh
cd ../superq
sh makesupq.sh
cd ../ripple
sh makeripp.sh
cd ../colormap
sh makemaps.sh
cd ../csg
sh makecsg.sh
cd ../texture
sh maketext.sh
cd ../misc
sh makemisc.sh
cd ../height
sh makeheig.sh
cd ../sweep
sh makeswep.sh
cd ../poly
sh makepoly.sh
cd ../gridded
sh makegrid.sh
cd ../blobs
sh makeblob.sh
cd ../uvmaps
sh makeuv.sh
cd ../person
polyray people.pi -o people.tga -M 512
cd ../function
sh makefunc.sh
cd ../stones
sh makeston.sh
cd ../animate
sh makeanim.sh
cd ..
