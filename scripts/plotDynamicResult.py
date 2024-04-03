import numpy as np
import matplotlib.pyplot as plt
import random
import xlwt
import pandas as pd
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

def main():
    time_list = []
    x_list = []
    val_list = []

    result_dir = r'..\tmpInputOutput\output.txt'
    with open(result_dir,'r') as txt:
        lines = txt.readlines()
        for line in lines:
            info = line.split()
            if info[0] == 'step':
                x_list_step = []
                val_list_step = []
                time_list.append(float(info[3]))
            if info[0][0] == 'x':
                x_list_step.append(float(info[1]))
                val_list_step.append(float(info[3]))
            if info[0][1:] == '21':
                x_list.append(x_list_step)
                val_list.append(val_list_step)

    for i in range(len(time_list)):
        save_dir = fr'..\figures\dynamic\dt2_{time_list[i]}.jpg'
        fig, ax = plt.subplots()
        ax.set_title(f"1D Wave Dynamic: Time {time_list[i]}s ")
        ax.set_ylim([-1.0,1.0])
        ax.set_xlim([0.0,20])
        ax.set_xlabel('x')
        ax.set_ylabel('u')
        ax.plot(x_list[i],val_list[i])
        #plt.show()
        plt.savefig(save_dir)


if __name__ == "__main__":
    main()