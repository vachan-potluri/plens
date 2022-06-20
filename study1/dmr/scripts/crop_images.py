# A simple script to crop images in "plots" directory. This enables uniform cropping for all
# figures. PIL library is used. For pixel coordinates, view the original figure in an image viewer
# that can display pixel location too (eg `geeqie`).

from PIL import Image

plot_dir = "/home/vachan/Documents/Work/plens/study1/dmr/plots" # loc for saving plots
for N in [1,2,3,5]:
    im = Image.open("{}/res_60_N{}_chandrashekhar.png".format(plot_dir, N))
    # remember image pixel coordinate axes start at the top left corner
    im_crop = im.crop((25, 186, 1420, 640)) # left, top, right, bottom
    im_crop.show()
    im_crop.save("{}/res_60_N{}_chandrashekhar_crop.png".format(plot_dir, N))