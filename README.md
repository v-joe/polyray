# polyray

Updated source code for Polyray raytracer by Alexander R. Enzmann that compiles under modern POSIX systems with X11 libraries (e.g. Linux, OSX, or Cygwin).

Original Polyray 1.8 source code and examples from https://web.archive.org/web/20030604183639/http://pages.infinit.net:80/gollum/polyray/poly1.htm.

## Build instructions

With pkg-config and X11 headers and libraries installed, running `make` is sufficient.<BR>
For OSX, make sure [XQuartz](https://github.com/XQuartz) and Xcode are installed and run `make -f makefile.osx`.

## Examples

The included examples (from the original Polyray distribution) can be regenerated in the `examples` directory with `sh makeall.sh`. `convert` from [ImageMagick](https://github.com/ImageMagick/ImageMagick) is required for generating the animation examples.

![image](https://raw.githubusercontent.com/v-joe/polyray/main/examples/animate/robarm/robarm.gif) &nbsp;&nbsp;&nbsp; ![image](https://raw.githubusercontent.com/v-joe/polyray/main/examples/animate/thresh/mvsph4.gif) &nbsp;&nbsp;&nbsp; ![image](https://raw.githubusercontent.com/v-joe/polyray/main/examples/animate/particle/part2.gif) &nbsp;&nbsp;&nbsp; ![image](https://raw.githubusercontent.com/v-joe/polyray/main/examples/animate/particle/part6.gif) &nbsp;&nbsp;&nbsp; ![image](https://raw.githubusercontent.com/v-joe/polyray/main/examples/animate/flying/flyplane.gif) &nbsp;&nbsp;&nbsp; ![image](https://raw.githubusercontent.com/v-joe/polyray/main/examples/animate/street/street.gif)
