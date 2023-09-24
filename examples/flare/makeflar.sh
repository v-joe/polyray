polyray lflare.pi  -o lflare.tga
polyray lflare2.pi -o lflare2.tga
polyray lflare3.pi -o lflare3.tga
polyray lflare4.pi -o lflare4.tga
polyray lflare5.pi -o lflare5.tga
convert out*.tga lflare5.gif
rm out*.tga
polyray flare0.pi  -o flare0.tga
polyray flare1.pi  -o flare1.tga
polyray flare2.pi
convert out*.tga flare2.gif
rm out*.tga
polyray flare3.pi
convert out*.tga flare3.gif
rm out*.tga
