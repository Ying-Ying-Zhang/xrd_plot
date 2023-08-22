import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import pymatgen as mg
from pymatgen.analysis.diffraction.xrd import XRDCalculator
import fnmatch
import re

def conv2cif(file_name, path='structures/',ftype='cif'):
    fname, ext = os.path.splitext(file_name)
    a = read(path_folder+file_name)
    out_file = fname + '.' + ftype
    write(path_folder+out_file, a, format=ftype)
    return out_file


def cal_xrd(file_name,two_theta_range=(0, 50)):
    fname, ext = os.path.splitext(file_name)
    if ext == 'gen':
        file_name = conv2cif(file_name,path=path)
    structure = mg.core.Structure.from_file(file_name)
    c = XRDCalculator(wavelength='CuKa1')
    xrd_data = XRDCalculator.get_pattern(c,structure,scaled=True,two_theta_range=two_theta_range)
    return xrd_data


def plot_xrd(x,y,name='generate_xrd.png', two_theta_range=(0,50)):
    fig, ax = plt.subplots(figsize=(6, 4),dpi=250)
    ax.plot(x, y)
    ax.set_xlabel("2θ(degrees)", fontsize=16)
    ax.xaxis.set_tick_params(labelsize=14, width=1.5)
    ax.yaxis.set_tick_params(labelsize=14, width=1.5)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.5)
    ax.set_xlim(two_theta_range)
    ax.set_ylabel("Intensity", fontsize=16)
    plt.tight_layout()
    plt.savefig(path_folder + '/' + name + '.png')


def spectrum(x,y,sigma,x_range):
    gE=[]
    for xi in x_range:
        tot=0
        for xj,os in zip(x,y):
            L=(FWHM/(2*np.pi))*(1/((xj-xi)**2+0.25*FWHM**2))
            G=a*np.exp(-(xj-xi)**2/(2*sigma**2))
            P=omega*G+(1-omega)*L
            tot+=os*P
            #tot+=os*np.exp(-((((xj-xi)/sigma)**2)))
        gE.append(tot)
    return gE


x_min = 0
x_max = 30
x_num = int((x_max - x_min) / 0.02) + 1
FWHM = 0.1
sigma = FWHM*0.42463  # 1/(2*sqrt(2*ln2))=0.42463
a = 1/(sigma*np.sqrt(2*np.pi))
omega = 0.001
print('sigma is:',sigma)
print('a is:',a)
x = np.linspace(x_min, x_max, num=x_num, endpoint=True)

path = 'C:/Downloads/'
#folder_list = ['IMDEA-COF-2_big-AA/', 'IMDEA-COF-2_small-AA/']
folder_list = os.listdir(path)

for i in folder_list:
    print(path+i)
    path_folder = path + i
    if os.path.exists(path_folder +'/generate_xrd.csv'):
        print('there is a generate_xrd.csv')
        f = pd.read_csv(path_folder +'/generate_xrd.csv')
        f.to_csv(path_folder +'/generate_xrd_bak.csv',index=False)
        if f.columns[-1] == 'intensity_ave':
            f = f.drop('intensity_ave', axis=1)
        print(re.split('_|.cif',f.columns[-1]))
        start_num = int(re.split('_|.cif',f.columns[-1])[1])+1
        intensity_sum = np.zeros([len(x)])
        for i in range(start_num):
            intensity_sum += f['50_'+str(i)+'.cif']
    else:
        print('no generate_xrd.csv')
        f = pd.DataFrame(x,columns=['twotheta'])
        f.to_csv(path_folder +'/generate_xrd.csv',index=False)
        intensity_sum = np.zeros([len(x)])
        start_num = 0
    files = fnmatch.filter(os.listdir(path_folder),'*.cif')
    num = len(files)
    print(start_num,num)
    print(files)
    for file in files:   # range(start_num, num):# files:
        #file = '50_'+str(j) + '.cif'  # j.lstrip()
        print('calculating', file)
        print(path_folder + '/' + file)
        xrd_data = cal_xrd(path_folder + '/'+ file,two_theta_range=(x_min,x_max))
        gxrd = spectrum(xrd_data.x, xrd_data.y, sigma, x)
        f[str(file)] = gxrd
        intensity_sum += np.array(gxrd)
        f.to_csv(path_folder + '/generate_xrd.csv',index=False)
        print(intensity_sum[0])
        plot_xrd(x, np.array(gxrd), name=str(file), two_theta_range=(x_min, x_max))
    print('calculating intensity_ave')

    #intensity_ave = intensity_sum/num
    #f['intensity_ave'] = intensity_ave
    f.to_csv(path_folder+'/generate_xrd.csv',index=False)
    #plot_xrd(x,intensity_ave,two_theta_range=(x_min,x_max))
#plt.show()

