[1mdiff --git a/Docs/Multiscat_HTML_Docs/.DS_Store b/Docs/Multiscat_HTML_Docs/.DS_Store[m
[1mdeleted file mode 100644[m
[1mindex 64948d6..0000000[m
Binary files a/Docs/Multiscat_HTML_Docs/.DS_Store and /dev/null differ
[1mdiff --git a/FourierLabels.in b/FourierLabels.in[m
[1mindex 26b59b4..9511227 100644[m
[1m--- a/FourierLabels.in[m
[1m+++ b/FourierLabels.in[m
[36m@@ -4094,3 +4094,4 @@[m
 29 31[m
 30 31[m
 31 31[m
[41m+[m
[1mdiff --git a/Makefile b/Makefile[m
[1mindex 3cfe769..2aaa4ed 100644[m
[1m--- a/Makefile[m
[1m+++ b/Makefile[m
[36m@@ -1,11 +1,12 @@[m
 #FFLAGS       = -real_size 128[m
 #FFLAGS      = -check_bounds -real_size 128[m
[31m-FFLAGS       = -O3 -mcmodel=large[m
[32m+[m[32mFFLAGS       = -O3 -mcmodel=large -fno-PIC#viv woz ere[m
 #FFLAGS       = -O[m
 #FFLAGS      = -check_bounds[m
 #Used for debugging[m
 #FFLAGS       = -O3 -mcmodel=large  -g -fbacktrace -Wall -fcheck=all[m
 [m
[32m+[m
 multiscat:		multiscat.o scatsub.o diagsub.o potsub.o[m
 			gfortran ${FFLAGS} -o multiscat multiscat.o scatsub.o diagsub.o potsub.o[m
 [m
[1mdiff --git a/Multiscat.conf b/Multiscat.conf[m
[1mindex 76d196e..3f8311a 100644[m
[1m--- a/Multiscat.conf[m
[1m+++ b/Multiscat.conf[m
[36m@@ -17,3 +17,4 @@[m [mscatCond.in	! The scattering conditions input file[m
 10001       !startindex[m
 10001       !endindex[m
 4       !helium mass[m
[41m+[m
[1mdiff --git a/diffrac10001.out b/diffrac10001.out[m
[1mindex ebcf672..00a2f42 100644[m
[1m--- a/diffrac10001.out[m
[1m+++ b/diffrac10001.out[m
[36m@@ -5,29 +5,29 @@[m
 # ------------------------------------[m
  # Beam energy           ei = 0.2000E+02[m
  # Polar angle        theta =   0.300000E+02  0.000000E+00[m
[31m-#     -4    -1          0.106161E-02[m
[31m-#     -4     0          0.139205E-02[m
[31m-#     -4     1          0.106161E-02[m
[31m-#     -3    -2          0.749058E-02[m
[31m-#     -3    -1          0.288746E-01[m
[31m-#     -3     0          0.206047E-01[m
[31m-#     -3     1          0.288746E-01[m
[31m-#     -3     2          0.749058E-02[m
[31m-#     -2    -2          0.473725E-01[m
[31m-#     -2    -1          0.831385E-01[m
[31m-#     -2     0          0.287208E-01[m
[31m-#     -2     1          0.831385E-01[m
[31m-#     -2     2          0.473725E-01[m
[31m-#     -1    -2          0.540619E-01[m
[31m-#     -1    -1          0.774560E-01[m
[31m-#     -1     0          0.111715E-01[m
[31m-#     -1     1          0.774560E-01[m
[31m-#     -1     2          0.540619E-01[m
[31m-#      0    -2          0.221086E-01[m
[31m-#      0    -1          0.206786E-01[m
[31m-#      0     0          0.245278E-01[m
[31m-#      0     1          0.206786E-01[m
[31m-#      0     2          0.221086E-01[m
[31m-#      1    -1          0.100728E+00[m
[31m-#      1     0          0.276419E-01[m
[31m-#      1     1          0.100728E+00[m
[32m+[m[32m#     -4    -1         -0.100000E+01[m
[32m+[m[32m#     -4     0         -0.100000E+01[m
[32m+[m[32m#     -4     1         -0.100000E+01[m
[32m+[m[32m#     -3    -2         -0.100000E+01[m
[32m+[m[32m#     -3    -1         -0.100000E+01[m
[32m+[m[32m#     -3     0         -0.100000E+01[m
[32m+[m[32m#     -3     1         -0.100000E+01[m
[32m+[m[32m#     -3     2         -0.100000E+01[m
[32m+[m[32m#     -2    -2         -0.100000E+01[m
[32m+[m[32m#     -2    -1         -0.100000E+01[m
[32m+[m[32m#     -2     0         -0.100000E+01[m
[32m+[m[32m#     -2     1         -0.100000E+01[m
[32m+[m[32m#     -2     2         -0.100000E+01[m
[32m+[m[32m#     -1    -2         -0.100000E+01[m
[32m+[m[32m#     -1    -1         -0.100000E+01[m
[32m+[m[32m#     -1     0         -0.100000E+01[m
[32m+[m[32m#     -1     1         -0.100000E+01[m
[32m+[m[32m#     -1     2         -0.100000E+01[m
[32m+[m[32m#      0    -2         -0.100000E+01[m
[32m+[m[32m#      0    -1         -0.100000E+01[m
[32m+[m[32m#      0     0         -0.100000E+01[m
[32m+[m[32m#      0     1         -0.100000E+01[m
[32m+[m[32m#      0     2         -0.100000E+01[m
[32m+[m[32m#      1    -1         -0.100000E+01[m
[32m+[m[32m#      1     0         -0.100000E+01[m
[32m+[m[32m#      1     1         -0.100000E+01[m
[1mdiff --git a/multiscat b/multiscat[m
[1mold mode 100644[m
[1mnew mode 100755[m
[1mindex 8e6bdb5..55a7827[m
Binary files a/multiscat and b/multiscat differ
[1mdiff --git a/multiscat_plot.py b/multiscat_plot.py[m
[1mindex f6e3918..51fd1cf 100644[m
[1m--- a/multiscat_plot.py[m
[1m+++ b/multiscat_plot.py[m
[36m@@ -1,3 +1,4 @@[m
[32m+[m[32m#%%[m
 # -*- coding: utf-8 -*-[m
 """[m
 Created on Tue Oct 10 16:42:27 2023[m
[36m@@ -10,7 +11,8 @@[m [mThis is a suggested method for plotting the output results of multiscat.[m
 import numpy as np[m
 import pandas as pd[m
 import seaborn as sns[m
[31m-[m
[32m+[m[32mimport datetime[m
[32m+[m[32mfrom matplotlib import pyplot as plt[m
 # Default theme[m
 sns.set_theme()[m
 [m
[36m@@ -28,3 +30,7 @@[m [max = sns.heatmap(d2, cmap='viridis', cbar_kws={'label' : '$P(n_1,n_2)$'})[m
 ax.set_aspect('equal')[m
 ax.set_xlabel('$n_1$')[m
 ax.set_ylabel('$n_2$')[m
[32m+[m
[32m+[m[32msavestr = "test/Figures/" + datetime.datetime.now().strftime('Diffraction_%Y-%m-%d_%H-%M') + ".png"[m
[32m+[m[32mplt.savefig(fname=savestr)[m
[32m+[m[32m# %%[m
[1mdiff --git a/pot10001.in b/pot10001.in[m
[1mindex 239ea71..6d2d945 100644[m
[1m--- a/pot10001.in[m
[1m+++ b/pot10001.in[m
[36m@@ -409603,3 +409603,4 @@[m [mDummy line5[m
 (+4.118687e-22, +2.055352e-22)[m
 (-1.545502e-23, +3.106222e-23)[m
 (+8.071476e-23, +1.333444e-22)[m
[41m+[m
