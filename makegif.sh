# NOTE: this overwrite the hires images with low res so be sure 
# to make a copy first
mogrify -resize 1920x1080 *png
convert -dispose previous -delay 70 *.png -delay 350 img0013.png tess_south.gif
