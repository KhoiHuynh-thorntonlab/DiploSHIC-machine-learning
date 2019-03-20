import subprocess
import tensorflow as tf
import os
from keras import backend as K
os.environ["CUDA_VISIBLE_DEVICES"]="0"


with tf.device('/device:GPU:0'):
    for i in range(0,30):
        for j in range(10001,10101):
            filename1= "/media/khoih/result/withcol/highmutratestat."+str(i)+"."+str(j)+".txt"
            filename2= "/media/khoih/result/result/highmutratestat."+str(i)+"."+str(j)+".result.txt"
            subprocess.call(["python3", "/media/khoih/diploSHIC/diploSHIC.py", "predict", "/media/khoih/model_11000th_2500_25000alpha/modelfile.txt.json", "/media/khoih/model_11000th_2500_25000alpha/modelfile.txt.weights.hdf5", filename1, filename2])

