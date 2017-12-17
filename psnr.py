import math
import numpy as np
import sys

refName = sys.argv[1]
imgName = sys.argv[2]

f1 = open(refName, 'rb')
f2 = open(imgName, 'rb')

ref = bytearray(f1.read())
img = bytearray(f2.read()) 

height = 512
width  = 512

mse = 0 
for y in range(height):
    for x in range(width):
        dif = ref[y*height + x] - img[y*height + x]
        mse = mse + math.pow(dif , 2) / (height * width)

psnr = 10 * math.log(pow(255, 2) / mse, 10)

print "PSNR:",psnr