#!/usr/bin/env /opt/local/bin/python3.6
# -*- coding: utf-8 -*-
#=================================
#   Module import
#=================================
try:
    import os
    import sys
    import glob
    import struct
    import math
    import re
    import subprocess
except ImportError:
    print("Module os,sys,glob,struct,math,re,subprocess are not found.")
    quit()

try:
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
except ImportError:
    print("Module matplotlib.pyplot,cm are not found.")
    quit()

try:
    import numpy as np
except ImportError:
    print("Module numpy is not found.")
    quit()

#================================
#   Class definition
#================================
class Plotter:
    def __init__(self,xmin,xmax,ymin,ymax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.x = []
        self.y = []
        self.z = []

    def read_file(self,filename,skip_freq):
        # Read a file
        print("reading {0} (skip_freq = {1})".format(filename,skip_freq))
        fp = open(filename,"r")
        data_num = 0
        for line in fp:
            items = line.split()
            try:
                x = float(items[0])
                y = float(items[1])
                z = float(items[2])
                if skip_freq == 0:
                    self.x.append(x)
                    self.y.append(y)
                    self.z.append(z)
                else:
                    data_num += 1
                    if (data_num == skip_freq):
                        self.x.append(x)
                        self.y.append(y)
                        self.z.append(z)
                        data_num = 0
                            
            except:
                print("cannot convert to FP values!")
                sys.exit()
        fp.close()
        print("{} doubles are read.".format(len(self.x)+len(self.y)+len(self.z)))

        # Convert to numpy-format
        self.x = np.array(self.x)
        self.y = np.array(self.y)
        self.z = np.array(self.z)
        print("converted to numpy-format.")

    def make_fig(self,basename,dpi,save_fmt='ps',make_xbb='no'):
        # Set file names
        ps_file   = basename + '.ps'
        eps_file  = basename + '.eps'
        pdf_file  = basename + '.pdf'
        png_file  = basename + '.png'
        
        # Make a figure
        font = {'family' : 'Verdana',
                'weight' : 'normal',
                'size'   : '14'}
        plt.rc('font',**font)
        plt.rc('text',usetex=True)
        #fig = plt.figure(1,figsize=(10,10))
        #ax = fig.add_subplot(111)
        #ax.set_aspect('equal')
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        
        # Plot column density of the reference models
        ax.set_xlabel(r'$x$',fontsize=20)
        ax.set_ylabel(r'$y$',fontsize=20)
        ax.tick_params(axis='x',labelsize=24)
        ax.tick_params(axis='y',labelsize=24)
        ax.set_xlim(self.xmin,self.xmax)
        ax.set_ylim(self.ymin,self.ymax)
        x_vals = np.resize(self.x,(np.size(self.x),1))
        y_vals = np.resize(self.y,(np.size(self.y),1))
        ax.scatter(x_vals, y_vals, s=0.16, c='k',edgecolors='none') 
        
        # Display figure 
        fig.tight_layout()

        # Save figure
        if (save_fmt == 'ps'):
            fig.savefig(ps_file)
            # Make eps,pdf,bb files
            cmd = 'ps2eps -B -l -g -a -f ' + ps_file
            subprocess.call(cmd,shell=True)
            cmd = 'eps2pdf -B -H -f ' + eps_file
            subprocess.call(cmd,shell=True)
            if (make_xbb == 'yes'):
                cmd = 'extractbb ' + pdf_file
                subprocess.call(cmd,shell=True)
        elif save_fmt == 'png':
            fig.savefig(png_file,dpi=dpi)
            if (make_xbb == 'yes'):
                cmd = 'extractbb ' + png_file
                subprocess.call(cmd,shell=True)
       
        # Close figure
        plt.close()

#================================
#   Main function
#================================
# Parameter settings
skip_freq = 0
dpi = 800
# Make figures of stellar distributions
#files = glob.glob("pos_star*txt")
files = glob.glob("pos_sph*txt")
for f in files:
    basename = os.path.splitext(os.path.basename(f))[0]
    xmin = - 0.1
    xmax =   0.1
    ymin = - 0.1
    ymax =   0.1
    P = Plotter(xmin,xmax,ymin,ymax)
    P.read_file(f,skip_freq)
    P.make_fig(basename, dpi, save_fmt='png')

