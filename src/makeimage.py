import sys
import os
from PIL import Image
import numpy as np
import csv
import math
import datetime



def dot3d(a,b):
	return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2])

def scale3d(a,scalar):
	return (a[0] * scalar, a[1] * scalar, a[2] * scalar)

def div3d(a,scalar):
	return (a[0] / scalar, a[1] / scalar, a[2] / scalar)

def round3d(a):
	return (int(round(a[0])), int(round(a[1])), int(round(a[2])))

DATA_DIR = 'data'
data_path = os.path.dirname(os.path.realpath(__file__)) + '/' + DATA_DIR


w = None
h = None

for csv_file_name in os.listdir(data_path):
	with open(data_path + "/" + csv_file_name, mode='r') as file:
		if not csv_file_name.endswith('.csv'):
			continue
		print( csv_file_name.split(':')[0:2])
		w, h = csv_file_name.split(':')[0:2]
		w = int(w)
		h = int(h)
		file.close()

assert w is not None

num_files = 0
pixels = []

for i in range(h *w):
	pixels.append([0,0,0])
[[0,0,0], [0,0,0], [0,0,0]] 



for csv_file_name in os.listdir(data_path):
	with open(data_path + "/" + csv_file_name, mode='r') as file:
		if not csv_file_name.endswith('.csv'):
			continue
		num_files += 1
		new_w, new_h, grey_conv, reinhard_tone_map_a = csv_file_name.split(':')[0:4]
		new_w = int(w)
		new_h = int(h)
		reinhard_tone_map_a = float(reinhard_tone_map_a)
		reinhard_tone_map_a /= 100
		assert new_h == h and new_w == w
		csv_reader = csv.reader(file, delimiter=',')
		for i, row in enumerate(csv_reader):
			if i == 0:
				continue
			pixels[i-1][0] += float(row[0])
			pixels[i-1][1] += float(row[1]) 
			pixels[i-1][2] += float(row[2])
		file.close()

gray_conv_coeffs = []
g1, g2, g3 = grey_conv.split('=')
gray_conv_coeffs.append(float(g1))
gray_conv_coeffs.append(float(g2))
gray_conv_coeffs.append(float(g3))

assert w == new_w
assert h == new_h


# get avg luminance
avg_lum = 0.0
pix_count = 0
for i in range(w):
	for j in range(h):
		pixels[pix_count] = div3d(pixels[pix_count], num_files)
		color = pixels[pix_count]
		curr_lum = dot3d(gray_conv_coeffs, color)
		avg_lum += math.log(sys.float_info.epsilon + curr_lum)
		pix_count += 1

avg_lum /= (h*w)
avg_lum = math.exp(avg_lum)

# get maximum scaled luminance
max_t_pix = 0.0
pix_count = 0
for i in range(w):
	for j in range(h):
		color = pixels[pix_count]
		curr_lum = dot3d(gray_conv_coeffs, color)
		t_pix = reinhard_tone_map_a * curr_lum / avg_lum
		if (t_pix > max_t_pix):
			max_t_pix = t_pix
		
		pix_count += 1


# apply operator
global_max_color = 0.0
pix_count = 0
for i in range(w):
	for j in range(h):
		color = pixels[pix_count]
		curr_lum = dot3d(gray_conv_coeffs, color)
		scaling = 1
		tone_map = 0
		if curr_lum > 0 + sys.float_info.epsilon:
			t_pix = reinhard_tone_map_a * curr_lum / avg_lum
			tone_map = t_pix * (1 + (t_pix / math.pow(max_t_pix, 2))) / (1 + t_pix)
			scaling = tone_map / curr_lum 

		color = scale3d(color, scaling)
		# get_max_color
		max_color =  max(color[0],color[1], color[2])      
		if max_color > global_max_color:
			global_max_color = max_color
		pixels[pix_count] = color
		pix_count += 1



image = Image.new("RGB", (w, h))
image_pix = image.load()

# normalize pixel values
pix_count = 0
for i in range(w):
	for j in reversed(range(h)):
		color = pixels[pix_count]
		color = div3d(color, global_max_color)
		print(color)
		pixels[pix_count] = color
		color = round3d(scale3d(color, 255))
		image_pix[i, j] = color
		pix_count += 1

now = datetime.datetime.now()
image.save("python" + str(now) + ".png")
print(num_files)
print(len(pixels))



	