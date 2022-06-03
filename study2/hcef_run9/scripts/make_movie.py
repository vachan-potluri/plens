# 03-Jun-2022
# Makes a movie out of streamline and pressure plots.
# Written for PMRF review presentation
# Follow the answer and comments here: https://stackoverflow.com/a/44948030/8028836

import cv2
import numpy as np

res_dir = "../trial3/temp/"
counters = [1736, 1757, 1777, 1797, 1820, 1830, 1840, 1860, 1880]
base_streamline_im_name = "T_surface_streamlines"
base_cp_im_name = "cp_comparison"
streamline_images = []
for c in counters:
    streamline_im_name = res_dir + base_streamline_im_name + "_{}.png".format(c)
    streamline_im = cv2.imread(streamline_im_name)
    cp_im_name = res_dir + base_cp_im_name + "_{}.png".format(c)
    cp_im = cv2.imread(cp_im_name)

    h1, w1, l1 = streamline_im.shape
    h2, w2, l2 = cp_im.shape
    h = max(h1, h2)

    # for vertical centering
    combined_im_mat = np.zeros((h, w1+w2, 3), np.uint8)
    combined_im_mat[int((h-h1)/2):h1+int((h-h1)/2), :w1, :3] = streamline_im
    combined_im_mat[int((h-h2)/2):h2+int((h-h2)/2), w1:w1+w2, :3] = cp_im
    
    size = (w1+w2, max(h1,h2))
    streamline_images.append(combined_im_mat)

fps = 1
fourcc = cv2.VideoWriter_fourcc('m', 'p', '4', 'v') # https://www.programcreek.com/python/example/72134/cv2.VideoWriter
out = cv2.VideoWriter(res_dir + "movie.mp4", fourcc, fps, size)
for im in streamline_images:
    out.write(im)
out.release()
