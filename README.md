# TESS Full Sky Mosaic
This repo contains the code used to create a full-sky mosaic of NASA's TESS 
full frame image data, like the gif shown below. Included in this is 
corner glow modeling to remove the increased background light observed in the
corners of the cameras and produce  an image with a more uniform background 
everywhere.

![TESS gif](tess_south_azeq_label.gif)


To recreate this gif, you can simply clone the repository and run
```bash
python full_sky.py
```

Required packages: `numpy, matplotlib, scipy, astropy, cartopy`.

Running the code will first download the 208 full frame images used in the 
mosaic (7 GB) into the data/ directory.

It will then create hi-res versions of the 13 frames of the gif (1.5 GB total) 
in the figs/gif_azeq_label/ directory. Warning: this will take several hours on 
a single machine and use all of your RAM.

If you have `imagemagick` installed, you can then recreate the gif. I recommend
copying the files into a lo-res directory first. Then you can navigate into the
appropriate directory and run the two commands in `makegif.sh` to downsize the
individual frames and create the gif.
