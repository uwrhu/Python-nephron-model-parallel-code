import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import os
import argparse

segment = ['PT','DL','mTAL','DCT','CNT','CCD','urine']

female_base = './Female_hum_normal'
female_normal_file = './ACE_Female_hum_normal'
female_nhe50_file = './ACE_Female_hum_pt' # nhe50/nkcc70
female_nhe80_file = './ACE_Female_hum_distal' # nhe80/nkcc100

male_base = './Male_hum_normal'
male_normal_file = './ACE_Male_hum_normal'
male_nhe50_file = './ACE_Male_hum_pt'
male_nhe80_file = './ACE_Male_hum_distal'

neph_weight = [0.85,(0.15)*0.4,(0.15)*0.3,(0.15)*0.15,(0.15)*0.1,(0.15)*0.05]

solute = ['Na','K','Cl','HCO3','H2CO3','CO2','HPO4','H2PO4','urea','NH3','NH4','H','HCO2','H2CO2','glu']
segment_early = ['pt','sdl','mtal','dct','cnt']
segment_late = ['ccd','imcd']

bar_width = 0.2
fig,axarr = plt.subplots(4,2)
fig.set_figheight(60)
fig.set_figwidth(40)
fig.subplots_adjust(hspace = 0.06)

volume_conversion = 1.44
solute_conversion = 1.44*10e-4

#=================================================================
# Na
#=================================================================
s = 'Na'

#==========================
# Baseline
#==========================

female_delivery_number_base = []
female_delivery_sup_base = []
female_delivery_jux1_base = []
female_delivery_jux2_base = []
female_delivery_jux3_base = []
female_delivery_jux4_base = []
female_delivery_jux5_base = []
for seg in segment_early:
    file_sup = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    female_delivery_number_base.append(0)
    female_delivery_sup_base.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    female_delivery_jux1_base.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    female_delivery_jux2_base.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    female_delivery_jux3_base.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    female_delivery_jux4_base.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    female_delivery_jux5_base.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        female_delivery_number_base.append(0)
        female_delivery_sup_base.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        female_delivery_jux1_base.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        female_delivery_jux2_base.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        female_delivery_jux3_base.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        female_delivery_jux4_base.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        female_delivery_jux5_base.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        female_delivery_number_base.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    

male_delivery_number_base = []
male_delivery_sup_base = []
male_delivery_jux1_base = []
male_delivery_jux2_base = []
male_delivery_jux3_base = []
male_delivery_jux4_base = []
male_delivery_jux5_base = []
for seg in segment_early:
    file_sup = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    male_delivery_number_base.append(0)
    male_delivery_sup_base.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    male_delivery_jux1_base.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    male_delivery_jux2_base.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    male_delivery_jux3_base.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    male_delivery_jux4_base.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    male_delivery_jux5_base.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        male_delivery_number_base.append(0)
        male_delivery_sup_base.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        male_delivery_jux1_base.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        male_delivery_jux2_base.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        male_delivery_jux3_base.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        male_delivery_jux4_base.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        male_delivery_jux5_base.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        male_delivery_number_base.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]

#==========================
female_delivery_number = []
female_delivery_sup = []
female_delivery_jux1 = []
female_delivery_jux2 = []
female_delivery_jux3 = []
female_delivery_jux4 = []
female_delivery_jux5 = []
for seg in segment_early:
    file_sup = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    female_delivery_number.append(0)
    female_delivery_sup.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    female_delivery_jux1.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    female_delivery_jux2.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    female_delivery_jux3.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    female_delivery_jux4.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    female_delivery_jux5.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        female_delivery_number.append(0)
        female_delivery_sup.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        female_delivery_jux1.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        female_delivery_jux2.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        female_delivery_jux3.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        female_delivery_jux4.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        female_delivery_jux5.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        female_delivery_number.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    

male_delivery_number = []
male_delivery_sup = []
male_delivery_jux1 = []
male_delivery_jux2 = []
male_delivery_jux3 = []
male_delivery_jux4 = []
male_delivery_jux5 = []
for seg in segment_early:
    file_sup = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    male_delivery_number.append(0)
    male_delivery_sup.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    male_delivery_jux1.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    male_delivery_jux2.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    male_delivery_jux3.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    male_delivery_jux4.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    male_delivery_jux5.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        male_delivery_number.append(0)
        male_delivery_sup.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        male_delivery_jux1.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        male_delivery_jux2.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        male_delivery_jux3.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        male_delivery_jux4.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        male_delivery_jux5.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        male_delivery_number.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    
    #===================================
    # NHE3 50% inhibited
    #===================================
female_delivery_number_nhe50 = []
female_delivery_sup_nhe50 = []
female_delivery_jux1_nhe50 = []
female_delivery_jux2_nhe50 = []
female_delivery_jux3_nhe50 = []
female_delivery_jux4_nhe50 = []
female_delivery_jux5_nhe50 = []
for seg in segment_early:
    file_sup = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    female_delivery_number_nhe50.append(0)
    female_delivery_sup_nhe50.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    female_delivery_jux1_nhe50.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    female_delivery_jux2_nhe50.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    female_delivery_jux3_nhe50.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    female_delivery_jux4_nhe50.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    female_delivery_jux5_nhe50.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        female_delivery_number_nhe50.append(0)
        female_delivery_sup_nhe50.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        female_delivery_jux1_nhe50.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        female_delivery_jux2_nhe50.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        female_delivery_jux3_nhe50.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        female_delivery_jux4_nhe50.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        female_delivery_jux5_nhe50.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        female_delivery_number_nhe50.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    

male_delivery_number_nhe50 = []
male_delivery_sup_nhe50 = []
male_delivery_jux1_nhe50 = []
male_delivery_jux2_nhe50 = []
male_delivery_jux3_nhe50 = []
male_delivery_jux4_nhe50 = []
male_delivery_jux5_nhe50 = []
for seg in segment_early:
    file_sup = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    male_delivery_number_nhe50.append(0)
    male_delivery_sup_nhe50.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    male_delivery_jux1_nhe50.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    male_delivery_jux2_nhe50.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    male_delivery_jux3_nhe50.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    male_delivery_jux4_nhe50.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    male_delivery_jux5_nhe50.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        male_delivery_number_nhe50.append(0)
        male_delivery_sup_nhe50.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        male_delivery_jux1_nhe50.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        male_delivery_jux2_nhe50.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        male_delivery_jux3_nhe50.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        male_delivery_jux4_nhe50.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        male_delivery_jux5_nhe50.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        male_delivery_number_nhe50.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    #===================================
    # NHE3 80% inhibited
    #===================================
female_delivery_number_nhe80 = []
female_delivery_sup_nhe80 = []
female_delivery_jux1_nhe80 = []
female_delivery_jux2_nhe80 = []
female_delivery_jux3_nhe80 = []
female_delivery_jux4_nhe80 = []
female_delivery_jux5_nhe80 = []
for seg in segment_early:
    file_sup = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    female_delivery_number_nhe80.append(0)
    female_delivery_sup_nhe80.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    female_delivery_jux1_nhe80.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    female_delivery_jux2_nhe80.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    female_delivery_jux3_nhe80.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    female_delivery_jux4_nhe80.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    female_delivery_jux5_nhe80.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        female_delivery_number_nhe80.append(0)
        female_delivery_sup_nhe80.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        female_delivery_jux1_nhe80.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        female_delivery_jux2_nhe80.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        female_delivery_jux3_nhe80.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        female_delivery_jux4_nhe80.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        female_delivery_jux5_nhe80.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        female_delivery_number_nhe80.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    

male_delivery_number_nhe80 = []
male_delivery_sup_nhe80 = []
male_delivery_jux1_nhe80 = []
male_delivery_jux2_nhe80 = []
male_delivery_jux3_nhe80 = []
male_delivery_jux4_nhe80 = []
male_delivery_jux5_nhe80 = []
for seg in segment_early:
    file_sup = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    male_delivery_number_nhe80.append(0)
    male_delivery_sup_nhe80.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    male_delivery_jux1_nhe80.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    male_delivery_jux2_nhe80.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    male_delivery_jux3_nhe80.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    male_delivery_jux4_nhe80.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    male_delivery_jux5_nhe80.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        male_delivery_number_nhe80.append(0)
        male_delivery_sup_nhe80.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        male_delivery_jux1_nhe80.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        male_delivery_jux2_nhe80.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        male_delivery_jux3_nhe80.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        male_delivery_jux4_nhe80.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        male_delivery_jux5_nhe80.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        male_delivery_number_nhe80.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]

male_sup_base=axarr[0,0].bar(np.arange(len(segment[:6]))-bar_width,male_delivery_sup_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='blue',label='Male baseline')
male_jux_base=axarr[0,0].bar(np.arange(len(segment[:6]))-bar_width,[male_delivery_jux1_base[i]+male_delivery_jux2_base[i]+male_delivery_jux3_base[i]+male_delivery_jux4_base[i]+male_delivery_jux5_base[i] for i in range(len(male_delivery_sup_base))],bar_width,bottom=male_delivery_sup_base,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_base=axarr[0,0].bar(np.arange(len(segment))-bar_width,male_delivery_number_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='blue')

male_sup=axarr[0,0].bar(np.arange(len(segment[:6])),male_delivery_sup,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='royalblue',label='Male full ACEi')
male_jux=axarr[0,0].bar(np.arange(len(segment[:6])),[male_delivery_jux1[i]+male_delivery_jux2[i]+male_delivery_jux3[i]+male_delivery_jux4[i]+male_delivery_jux5[i] for i in range(len(male_delivery_sup))],bar_width,bottom=male_delivery_sup,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later=axarr[0,0].bar(np.arange(len(segment)),male_delivery_number,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='royalblue')

male_sup_nhe50=axarr[0,0].bar(np.arange(len(segment[:6]))+bar_width,male_delivery_sup_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='deepskyblue',label='Male proximal ACEi')
male_jux_nhe50=axarr[0,0].bar(np.arange(len(segment[:6]))+bar_width,[male_delivery_jux1_nhe50[i]+male_delivery_jux2_nhe50[i]+male_delivery_jux3_nhe50[i]+male_delivery_jux4_nhe50[i]+male_delivery_jux5_nhe50[i] for i in range(len(male_delivery_sup_nhe50))],bar_width,bottom=male_delivery_sup_nhe50,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_nhe50=axarr[0,0].bar(np.arange(len(segment))+bar_width,male_delivery_number_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='deepskyblue')

male_sup_nhe80=axarr[0,0].bar(np.arange(len(segment[:6]))+2*bar_width,male_delivery_sup_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='paleturquoise',label='Male distal ACEi')
male_jux_nhe80=axarr[0,0].bar(np.arange(len(segment[:6]))+2*bar_width,[male_delivery_jux1_nhe80[i]+male_delivery_jux2_nhe80[i]+male_delivery_jux3_nhe80[i]+male_delivery_jux4_nhe80[i]+male_delivery_jux5_nhe80[i] for i in range(len(male_delivery_sup_nhe80))],bar_width,bottom=male_delivery_sup_nhe80,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_nhe80=axarr[0,0].bar(np.arange(len(segment))+2*bar_width,male_delivery_number_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='paleturquoise')

female_sup_base=axarr[0,1].bar(np.arange(len(segment[:6]))-bar_width,female_delivery_sup_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='red',label='Female baseline')
female_jux_base=axarr[0,1].bar(np.arange(len(segment[:6]))-bar_width,[female_delivery_jux1_base[i]+female_delivery_jux2_base[i]+female_delivery_jux3_base[i]+female_delivery_jux4_base[i]+female_delivery_jux5_base[i] for i in range(len(female_delivery_sup_base))],bar_width,bottom=female_delivery_sup_base,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
female_later_base=axarr[0,1].bar(np.arange(len(segment))-bar_width,female_delivery_number_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='red')

Female_sup=axarr[0,1].bar(np.arange(len(segment[:6])),female_delivery_sup,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='magenta',label='Female full ACEi')
Female_jux=axarr[0,1].bar(np.arange(len(segment[:6])),[female_delivery_jux1[i]+female_delivery_jux2[i]+female_delivery_jux3[i]+female_delivery_jux4[i]+female_delivery_jux5[i] for i in range(len(female_delivery_sup))],bar_width,bottom=female_delivery_sup,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later=axarr[0,1].bar(np.arange(len(segment)),female_delivery_number,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='magenta')

Female_sup_nhe50=axarr[0,1].bar(np.arange(len(segment[:6]))+1*bar_width,female_delivery_sup_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='hotpink',label='Female proximal ACEi')
Female_jux_nhe50=axarr[0,1].bar(np.arange(len(segment[:6]))+1*bar_width,[female_delivery_jux1_nhe50[i]+female_delivery_jux2_nhe50[i]+female_delivery_jux3_nhe50[i]+female_delivery_jux4_nhe50[i]+female_delivery_jux5_nhe50[i] for i in range(len(female_delivery_sup_nhe50))],bar_width,bottom=female_delivery_sup_nhe50,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later_nhe50=axarr[0,1].bar(np.arange(len(segment))+1*bar_width,female_delivery_number_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='hotpink')

Female_sup_nhe80=axarr[0,1].bar(np.arange(len(segment[:6]))+2*bar_width,female_delivery_sup_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='pink',label='Female distal ACEi')
Female_jux_nhe80=axarr[0,1].bar(np.arange(len(segment[:6]))+2*bar_width,[female_delivery_jux1_nhe80[i]+female_delivery_jux2_nhe80[i]+female_delivery_jux3_nhe80[i]+female_delivery_jux4_nhe80[i]+female_delivery_jux5_nhe80[i] for i in range(len(female_delivery_sup_nhe80))],bar_width,bottom=female_delivery_sup_nhe80,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later_nhe80=axarr[0,1].bar(np.arange(len(segment))+2*bar_width,female_delivery_number_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='pink')

axarr[0,0].set_xticks(np.arange(len(segment))+0.5*bar_width)
axarr[0,0].set_xticklabels(segment,fontsize=40)
axarr[0,0].tick_params(axis='both',labelsize=40)
#axarr[0,0].set_ylim(0,185)
axarr[0,0].set_title('Male',fontsize = 50)
#ax.set_xlabel('Segment',fontsize=20)
axarr[0,0].set_ylabel('Na$^+$ delivery (mol/Day)',fontsize=40)
axarr[0,0].legend(fontsize=30,markerscale=30)

axarr[0,1].set_xticks(np.arange(len(segment))+0.5*bar_width)
axarr[0,1].set_xticklabels(segment,fontsize=40)
axarr[0,1].tick_params(axis='both',labelsize=40)
#axarr[0,1].set_ylim(0,185)
axarr[0,1].set_title('Female',fontsize = 50)
#ax.set_xlabel('Segment',fontsize=20)
axarr[0,1].set_ylabel('Na$^+$ delivery (mol/Day)',fontsize=40)
axarr[0,1].legend(fontsize=30,markerscale=30)

bar_width_ins = bar_width
axins = inset_axes(axarr[0,0],width=3.5,height=3.5,loc=7)

male_sup_base_inset=axins.bar(np.arange(len(segment[:6]))-bar_width,male_delivery_sup_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='blue',label='Men, Baseline')
male_jux_base_inset=axins.bar(np.arange(len(segment[:6]))-bar_width,[male_delivery_jux1_base[i]+male_delivery_jux2_base[i]+male_delivery_jux3_base[i]+male_delivery_jux4_base[i]+male_delivery_jux5_base[i] for i in range(len(male_delivery_sup_base))],bar_width,bottom=male_delivery_sup_base,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_base_inset=axins.bar(np.arange(len(segment))-bar_width,male_delivery_number_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='blue')

male_sup_inset=axins.bar(np.arange(len(segment[:6])),male_delivery_sup,bar_width_ins,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='royalblue',label='Men, ACEi')
male_jux_inset=axins.bar(np.arange(len(segment[:6])),[male_delivery_jux1[i]+male_delivery_jux2[i]+male_delivery_jux3[i]+male_delivery_jux4[i]+male_delivery_jux5[i] for i in range(len(male_delivery_sup))],bar_width_ins,bottom=male_delivery_sup,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_inset=axins.bar(np.arange(len(segment)),male_delivery_number,bar_width_ins,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='royalblue')

male_sup_nhe50_inset=axins.bar(np.arange(len(segment[:6]))+bar_width,male_delivery_sup_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='deepskyblue',label='Men, ACEi PT')
male_jux_nhe50_inset=axins.bar(np.arange(len(segment[:6]))+bar_width,[male_delivery_jux1_nhe50[i]+male_delivery_jux2_nhe50[i]+male_delivery_jux3_nhe50[i]+male_delivery_jux4_nhe50[i]+male_delivery_jux5_nhe50[i] for i in range(len(male_delivery_sup_nhe50))],bar_width,bottom=male_delivery_sup_nhe50,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_nhe50_inset=axins.bar(np.arange(len(segment))+bar_width,male_delivery_number_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='deepskyblue')

male_sup_nhe80_inset=axins.bar(np.arange(len(segment[:6]))+2*bar_width,male_delivery_sup_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='paleturquoise',label='Men, ACEi distal')
male_jux_nhe80_inset=axins.bar(np.arange(len(segment[:6]))+2*bar_width,[male_delivery_jux1_nhe80[i]+male_delivery_jux2_nhe80[i]+male_delivery_jux3_nhe80[i]+male_delivery_jux4_nhe80[i]+male_delivery_jux5_nhe80[i] for i in range(len(male_delivery_sup_nhe80))],bar_width,bottom=male_delivery_sup_nhe80,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_nhe80_inset=axins.bar(np.arange(len(segment))+2*bar_width,male_delivery_number_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='paleturquoise')

axins.set_xticks(np.arange(len(segment))+0.5*bar_width_ins)
axins.set_xticklabels(segment,fontsize=40)
axins.set_xlim(5-2.5*bar_width_ins,6+3*bar_width_ins)
axins.set_ylim(0,1)
axins.tick_params(axis='both',labelsize=40)

bar_width_ins = bar_width
axins = inset_axes(axarr[0,1],width=3.5,height=3.5,loc=7)

female_sup_base_inset=axins.bar(np.arange(len(segment[:6]))-bar_width,female_delivery_sup_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='red',label='Women, Baseline')
female_jux_base_inset=axins.bar(np.arange(len(segment[:6]))-bar_width,[female_delivery_jux1_base[i]+female_delivery_jux2_base[i]+female_delivery_jux3_base[i]+female_delivery_jux4_base[i]+female_delivery_jux5_base[i] for i in range(len(female_delivery_sup_base))],bar_width,bottom=female_delivery_sup_base,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
female_later_base_inset=axins.bar(np.arange(len(segment))-bar_width,female_delivery_number_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='red')

Female_sup_inset=axins.bar(np.arange(len(segment[:6])),female_delivery_sup,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='magenta',label='Baseline')
Female_jux_inset=axins.bar(np.arange(len(segment[:6])),[female_delivery_jux1[i]+female_delivery_jux2[i]+female_delivery_jux3[i]+female_delivery_jux4[i]+female_delivery_jux5[i] for i in range(len(female_delivery_sup))],bar_width,bottom=female_delivery_sup,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later_inset=axins.bar(np.arange(len(segment)),female_delivery_number,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='magenta')

Female_sup_nhe50_inset=axins.bar(np.arange(len(segment[:6]))+1*bar_width,female_delivery_sup_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='hotpink',label='50% NHE3 inhib')
Female_jux_nhe50_inset=axins.bar(np.arange(len(segment[:6]))+1*bar_width,[female_delivery_jux1_nhe50[i]+female_delivery_jux2_nhe50[i]+female_delivery_jux3_nhe50[i]+female_delivery_jux4_nhe50[i]+female_delivery_jux5_nhe50[i] for i in range(len(female_delivery_sup_nhe50))],bar_width,bottom=female_delivery_sup_nhe50,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later_nhe50_inset=axins.bar(np.arange(len(segment))+1*bar_width,female_delivery_number_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='hotpink')

Female_sup_nhe80_inset=axins.bar(np.arange(len(segment[:6]))+2*bar_width,female_delivery_sup_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='pink',label='80% NHE3 inhib')
Female_jux_nhe80_inset=axins.bar(np.arange(len(segment[:6]))+2*bar_width,[female_delivery_jux1_nhe80[i]+female_delivery_jux2_nhe80[i]+female_delivery_jux3_nhe80[i]+female_delivery_jux4_nhe80[i]+female_delivery_jux5_nhe80[i] for i in range(len(female_delivery_sup_nhe80))],bar_width,bottom=female_delivery_sup_nhe80,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later_nhe80_inset=axins.bar(np.arange(len(segment))+2*bar_width,female_delivery_number_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='pink')

axins.set_xticks(np.arange(len(segment))+0.5*bar_width_ins)
axins.set_xticklabels(segment,fontsize=40)
axins.set_xlim(5-2.5*bar_width_ins,6+3*bar_width_ins)
axins.set_ylim(0,1)
axins.tick_params(axis='both',labelsize=40)



#=================================================================
# K
#=================================================================

s = 'K'
#==========================
# Baseline
#==========================

female_delivery_number_base = []
female_delivery_sup_base = []
female_delivery_jux1_base = []
female_delivery_jux2_base = []
female_delivery_jux3_base = []
female_delivery_jux4_base = []
female_delivery_jux5_base = []
for seg in segment_early:
    file_sup = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    female_delivery_number_base.append(0)
    female_delivery_sup_base.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    female_delivery_jux1_base.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    female_delivery_jux2_base.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    female_delivery_jux3_base.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    female_delivery_jux4_base.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    female_delivery_jux5_base.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        female_delivery_number_base.append(0)
        female_delivery_sup_base.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        female_delivery_jux1_base.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        female_delivery_jux2_base.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        female_delivery_jux3_base.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        female_delivery_jux4_base.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        female_delivery_jux5_base.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        female_delivery_number_base.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    

male_delivery_number_base = []
male_delivery_sup_base = []
male_delivery_jux1_base = []
male_delivery_jux2_base = []
male_delivery_jux3_base = []
male_delivery_jux4_base = []
male_delivery_jux5_base = []
for seg in segment_early:
    file_sup = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    male_delivery_number_base.append(0)
    male_delivery_sup_base.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    male_delivery_jux1_base.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    male_delivery_jux2_base.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    male_delivery_jux3_base.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    male_delivery_jux4_base.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    male_delivery_jux5_base.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        male_delivery_number_base.append(0)
        male_delivery_sup_base.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        male_delivery_jux1_base.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        male_delivery_jux2_base.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        male_delivery_jux3_base.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        male_delivery_jux4_base.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        male_delivery_jux5_base.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        male_delivery_number_base.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]

#==========================
female_delivery_number = []
female_delivery_sup = []
female_delivery_jux1 = []
female_delivery_jux2 = []
female_delivery_jux3 = []
female_delivery_jux4 = []
female_delivery_jux5 = []
for seg in segment_early:
    file_sup = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    female_delivery_number.append(0)
    female_delivery_sup.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    female_delivery_jux1.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    female_delivery_jux2.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    female_delivery_jux3.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    female_delivery_jux4.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    female_delivery_jux5.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        female_delivery_number.append(0)
        female_delivery_sup.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        female_delivery_jux1.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        female_delivery_jux2.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        female_delivery_jux3.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        female_delivery_jux4.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        female_delivery_jux5.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        female_delivery_number.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    

male_delivery_number = []
male_delivery_sup = []
male_delivery_jux1 = []
male_delivery_jux2 = []
male_delivery_jux3 = []
male_delivery_jux4 = []
male_delivery_jux5 = []
for seg in segment_early:
    file_sup = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    male_delivery_number.append(0)
    male_delivery_sup.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    male_delivery_jux1.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    male_delivery_jux2.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    male_delivery_jux3.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    male_delivery_jux4.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    male_delivery_jux5.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        male_delivery_number.append(0)
        male_delivery_sup.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        male_delivery_jux1.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        male_delivery_jux2.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        male_delivery_jux3.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        male_delivery_jux4.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        male_delivery_jux5.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        male_delivery_number.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    
    #===================================
    # NHE3 50% inhibited
    #===================================
female_delivery_number_nhe50 = []
female_delivery_sup_nhe50 = []
female_delivery_jux1_nhe50 = []
female_delivery_jux2_nhe50 = []
female_delivery_jux3_nhe50 = []
female_delivery_jux4_nhe50 = []
female_delivery_jux5_nhe50 = []
for seg in segment_early:
    file_sup = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    female_delivery_number_nhe50.append(0)
    female_delivery_sup_nhe50.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    female_delivery_jux1_nhe50.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    female_delivery_jux2_nhe50.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    female_delivery_jux3_nhe50.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    female_delivery_jux4_nhe50.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    female_delivery_jux5_nhe50.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        female_delivery_number_nhe50.append(0)
        female_delivery_sup_nhe50.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        female_delivery_jux1_nhe50.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        female_delivery_jux2_nhe50.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        female_delivery_jux3_nhe50.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        female_delivery_jux4_nhe50.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        female_delivery_jux5_nhe50.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        female_delivery_number_nhe50.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    

male_delivery_number_nhe50 = []
male_delivery_sup_nhe50 = []
male_delivery_jux1_nhe50 = []
male_delivery_jux2_nhe50 = []
male_delivery_jux3_nhe50 = []
male_delivery_jux4_nhe50 = []
male_delivery_jux5_nhe50 = []
for seg in segment_early:
    file_sup = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    male_delivery_number_nhe50.append(0)
    male_delivery_sup_nhe50.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    male_delivery_jux1_nhe50.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    male_delivery_jux2_nhe50.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    male_delivery_jux3_nhe50.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    male_delivery_jux4_nhe50.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    male_delivery_jux5_nhe50.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        male_delivery_number_nhe50.append(0)
        male_delivery_sup_nhe50.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        male_delivery_jux1_nhe50.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        male_delivery_jux2_nhe50.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        male_delivery_jux3_nhe50.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        male_delivery_jux4_nhe50.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        male_delivery_jux5_nhe50.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        male_delivery_number_nhe50.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    #===================================
    # NHE3 80% inhibited
    #===================================
female_delivery_number_nhe80 = []
female_delivery_sup_nhe80 = []
female_delivery_jux1_nhe80 = []
female_delivery_jux2_nhe80 = []
female_delivery_jux3_nhe80 = []
female_delivery_jux4_nhe80 = []
female_delivery_jux5_nhe80 = []
for seg in segment_early:
    file_sup = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    female_delivery_number_nhe80.append(0)
    female_delivery_sup_nhe80.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    female_delivery_jux1_nhe80.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    female_delivery_jux2_nhe80.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    female_delivery_jux3_nhe80.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    female_delivery_jux4_nhe80.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    female_delivery_jux5_nhe80.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        female_delivery_number_nhe80.append(0)
        female_delivery_sup_nhe80.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        female_delivery_jux1_nhe80.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        female_delivery_jux2_nhe80.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        female_delivery_jux3_nhe80.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        female_delivery_jux4_nhe80.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        female_delivery_jux5_nhe80.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        female_delivery_number_nhe80.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    

male_delivery_number_nhe80 = []
male_delivery_sup_nhe80 = []
male_delivery_jux1_nhe80 = []
male_delivery_jux2_nhe80 = []
male_delivery_jux3_nhe80 = []
male_delivery_jux4_nhe80 = []
male_delivery_jux5_nhe80 = []
for seg in segment_early:
    file_sup = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    male_delivery_number_nhe80.append(0)
    male_delivery_sup_nhe80.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    male_delivery_jux1_nhe80.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    male_delivery_jux2_nhe80.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    male_delivery_jux3_nhe80.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    male_delivery_jux4_nhe80.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    male_delivery_jux5_nhe80.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        male_delivery_number_nhe80.append(0)
        male_delivery_sup_nhe80.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        male_delivery_jux1_nhe80.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        male_delivery_jux2_nhe80.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        male_delivery_jux3_nhe80.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        male_delivery_jux4_nhe80.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        male_delivery_jux5_nhe80.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        male_delivery_number_nhe80.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]

male_sup_base=axarr[1,0].bar(np.arange(len(segment[:6]))-bar_width,male_delivery_sup_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='blue',label='Men, Baseline')
male_jux_base=axarr[1,0].bar(np.arange(len(segment[:6]))-bar_width,[male_delivery_jux1_base[i]+male_delivery_jux2_base[i]+male_delivery_jux3_base[i]+male_delivery_jux4_base[i]+male_delivery_jux5_base[i] for i in range(len(male_delivery_sup_base))],bar_width,bottom=male_delivery_sup_base,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_base=axarr[1,0].bar(np.arange(len(segment))-bar_width,male_delivery_number_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='blue')

male_sup=axarr[1,0].bar(np.arange(len(segment[:6])),male_delivery_sup,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='royalblue',label='Men, ACEi')
male_jux=axarr[1,0].bar(np.arange(len(segment[:6])),[male_delivery_jux1[i]+male_delivery_jux2[i]+male_delivery_jux3[i]+male_delivery_jux4[i]+male_delivery_jux5[i] for i in range(len(male_delivery_sup))],bar_width,bottom=male_delivery_sup,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later=axarr[1,0].bar(np.arange(len(segment)),male_delivery_number,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='royalblue')

male_sup_nhe50=axarr[1,0].bar(np.arange(len(segment[:6]))+bar_width,male_delivery_sup_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='deepskyblue',label='Men, NHE3 PT alter')
male_jux_nhe50=axarr[1,0].bar(np.arange(len(segment[:6]))+bar_width,[male_delivery_jux1_nhe50[i]+male_delivery_jux2_nhe50[i]+male_delivery_jux3_nhe50[i]+male_delivery_jux4_nhe50[i]+male_delivery_jux5_nhe50[i] for i in range(len(male_delivery_sup_nhe50))],bar_width,bottom=male_delivery_sup_nhe50,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_nhe50=axarr[1,0].bar(np.arange(len(segment))+bar_width,male_delivery_number_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='deepskyblue')

male_sup_nhe80=axarr[1,0].bar(np.arange(len(segment[:6]))+2*bar_width,male_delivery_sup_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='paleturquoise',label='Men, distal alter')
male_jux_nhe80=axarr[1,0].bar(np.arange(len(segment[:6]))+2*bar_width,[male_delivery_jux1_nhe80[i]+male_delivery_jux2_nhe80[i]+male_delivery_jux3_nhe80[i]+male_delivery_jux4_nhe80[i]+male_delivery_jux5_nhe80[i] for i in range(len(male_delivery_sup_nhe80))],bar_width,bottom=male_delivery_sup_nhe80,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_nhe80=axarr[1,0].bar(np.arange(len(segment))+2*bar_width,male_delivery_number_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='paleturquoise')

female_sup_base=axarr[1,1].bar(np.arange(len(segment[:6]))-bar_width,female_delivery_sup_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='red',label='Women, Baseline')
female_jux_base=axarr[1,1].bar(np.arange(len(segment[:6]))-bar_width,[female_delivery_jux1_base[i]+female_delivery_jux2_base[i]+female_delivery_jux3_base[i]+female_delivery_jux4_base[i]+female_delivery_jux5_base[i] for i in range(len(female_delivery_sup_base))],bar_width,bottom=female_delivery_sup_base,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
female_later_base=axarr[1,1].bar(np.arange(len(segment))-bar_width,female_delivery_number_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='red')

Female_sup=axarr[1,1].bar(np.arange(len(segment[:6])),female_delivery_sup,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='magenta',label='Women, ACEi')
Female_jux=axarr[1,1].bar(np.arange(len(segment[:6])),[female_delivery_jux1[i]+female_delivery_jux2[i]+female_delivery_jux3[i]+female_delivery_jux4[i]+female_delivery_jux5[i] for i in range(len(female_delivery_sup))],bar_width,bottom=female_delivery_sup,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later=axarr[1,1].bar(np.arange(len(segment)),female_delivery_number,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='magenta')

Female_sup_nhe50=axarr[1,1].bar(np.arange(len(segment[:6]))+1*bar_width,female_delivery_sup_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='hotpink',label='Women, NHE3 PT alter')
Female_jux_nhe50=axarr[1,1].bar(np.arange(len(segment[:6]))+1*bar_width,[female_delivery_jux1_nhe50[i]+female_delivery_jux2_nhe50[i]+female_delivery_jux3_nhe50[i]+female_delivery_jux4_nhe50[i]+female_delivery_jux5_nhe50[i] for i in range(len(female_delivery_sup_nhe50))],bar_width,bottom=female_delivery_sup_nhe50,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later_nhe50=axarr[1,1].bar(np.arange(len(segment))+1*bar_width,female_delivery_number_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='hotpink')

Female_sup_nhe80=axarr[1,1].bar(np.arange(len(segment[:6]))+2*bar_width,female_delivery_sup_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='pink',label='Women, distal alter')
Female_jux_nhe80=axarr[1,1].bar(np.arange(len(segment[:6]))+2*bar_width,[female_delivery_jux1_nhe80[i]+female_delivery_jux2_nhe80[i]+female_delivery_jux3_nhe80[i]+female_delivery_jux4_nhe80[i]+female_delivery_jux5_nhe80[i] for i in range(len(female_delivery_sup_nhe80))],bar_width,bottom=female_delivery_sup_nhe80,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later_nhe80=axarr[1,1].bar(np.arange(len(segment))+2*bar_width,female_delivery_number_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='pink')

axarr[1,0].set_xticks(np.arange(len(segment))+1*bar_width)
axarr[1,0].set_xticklabels(segment,fontsize=40)
axarr[1,0].tick_params(axis='both',labelsize=40)
#axarr[1,0].set_ylim(0,7)
#ax.set_xlabel('Segment',fontsize=20)
axarr[1,0].set_ylabel('K$^+$ delivery (mol/Day)',fontsize=40)
#axarr[1,0].legend(fontsize=30,markerscale=30)
bar_width_ins = bar_width
#axins = inset_axes(ax,width=2,height=2,loc=7)
#Male_inset=axins.bar(np.arange(len(segment)),male_list,bar_width_ins,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='black',label='Male')
#Female_inset=axins.bar(np.arange(len(segment))+bar_width_ins,female_list,bar_width_ins,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white',label='Female')
#axins.set_xticks(np.arange(len(segment))+0.5*bar_width_ins)
#axins.set_xticklabels(segment,fontsize=20)
#axins.set_xlim(5-1.5*bar_width_ins,6+2*bar_width)
#axins.set_ylim(0,2)
#axins.tick_params(axis='both',labelsize=20)
#autolabel(Male,'left')
#autolabel(Female,'right')
axarr[1,1].set_xticks(np.arange(len(segment))+1*bar_width)
axarr[1,1].set_xticklabels(segment,fontsize=40)
axarr[1,1].tick_params(axis='both',labelsize=40)
#axarr[1,1].set_ylim(0,7)
#ax.set_xlabel('Segment',fontsize=20)
axarr[1,1].set_ylabel('K$^+$ delivery (mol/Day)',fontsize=40)
#axarr[1,1].legend(fontsize=30,markerscale=30)
bar_width_ins = bar_width

#=================================================================
# Cl
#=================================================================

s = 'Cl'
#==========================
# Baseline
#==========================

female_delivery_number_base = []
female_delivery_sup_base = []
female_delivery_jux1_base = []
female_delivery_jux2_base = []
female_delivery_jux3_base = []
female_delivery_jux4_base = []
female_delivery_jux5_base = []
for seg in segment_early:
    file_sup = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    female_delivery_number_base.append(0)
    female_delivery_sup_base.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    female_delivery_jux1_base.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    female_delivery_jux2_base.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    female_delivery_jux3_base.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    female_delivery_jux4_base.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    female_delivery_jux5_base.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        female_delivery_number_base.append(0)
        female_delivery_sup_base.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        female_delivery_jux1_base.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        female_delivery_jux2_base.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        female_delivery_jux3_base.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        female_delivery_jux4_base.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        female_delivery_jux5_base.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(female_base+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        female_delivery_number_base.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    

male_delivery_number_base = []
male_delivery_sup_base = []
male_delivery_jux1_base = []
male_delivery_jux2_base = []
male_delivery_jux3_base = []
male_delivery_jux4_base = []
male_delivery_jux5_base = []
for seg in segment_early:
    file_sup = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    male_delivery_number_base.append(0)
    male_delivery_sup_base.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    male_delivery_jux1_base.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    male_delivery_jux2_base.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    male_delivery_jux3_base.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    male_delivery_jux4_base.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    male_delivery_jux5_base.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        male_delivery_number_base.append(0)
        male_delivery_sup_base.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        male_delivery_jux1_base.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        male_delivery_jux2_base.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        male_delivery_jux3_base.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        male_delivery_jux4_base.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        male_delivery_jux5_base.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(male_base+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        male_delivery_number_base.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]

#==========================
female_delivery_number = []
female_delivery_sup = []
female_delivery_jux1 = []
female_delivery_jux2 = []
female_delivery_jux3 = []
female_delivery_jux4 = []
female_delivery_jux5 = []
for seg in segment_early:
    file_sup = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    female_delivery_number.append(0)
    female_delivery_sup.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    female_delivery_jux1.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    female_delivery_jux2.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    female_delivery_jux3.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    female_delivery_jux4.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    female_delivery_jux5.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        female_delivery_number.append(0)
        female_delivery_sup.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        female_delivery_jux1.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        female_delivery_jux2.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        female_delivery_jux3.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        female_delivery_jux4.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        female_delivery_jux5.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(female_normal_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        female_delivery_number.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    

male_delivery_number = []
male_delivery_sup = []
male_delivery_jux1 = []
male_delivery_jux2 = []
male_delivery_jux3 = []
male_delivery_jux4 = []
male_delivery_jux5 = []
for seg in segment_early:
    file_sup = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    male_delivery_number.append(0)
    male_delivery_sup.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    male_delivery_jux1.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    male_delivery_jux2.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    male_delivery_jux3.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    male_delivery_jux4.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    male_delivery_jux5.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        male_delivery_number.append(0)
        male_delivery_sup.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        male_delivery_jux1.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        male_delivery_jux2.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        male_delivery_jux3.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        male_delivery_jux4.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        male_delivery_jux5.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(male_normal_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        male_delivery_number.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    
    #===================================
    # NHE3 50% inhibited
    #===================================
female_delivery_number_nhe50 = []
female_delivery_sup_nhe50 = []
female_delivery_jux1_nhe50 = []
female_delivery_jux2_nhe50 = []
female_delivery_jux3_nhe50 = []
female_delivery_jux4_nhe50 = []
female_delivery_jux5_nhe50 = []
for seg in segment_early:
    file_sup = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    female_delivery_number_nhe50.append(0)
    female_delivery_sup_nhe50.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    female_delivery_jux1_nhe50.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    female_delivery_jux2_nhe50.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    female_delivery_jux3_nhe50.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    female_delivery_jux4_nhe50.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    female_delivery_jux5_nhe50.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        female_delivery_number_nhe50.append(0)
        female_delivery_sup_nhe50.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        female_delivery_jux1_nhe50.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        female_delivery_jux2_nhe50.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        female_delivery_jux3_nhe50.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        female_delivery_jux4_nhe50.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        female_delivery_jux5_nhe50.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(female_nhe50_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        female_delivery_number_nhe50.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    

male_delivery_number_nhe50 = []
male_delivery_sup_nhe50 = []
male_delivery_jux1_nhe50 = []
male_delivery_jux2_nhe50 = []
male_delivery_jux3_nhe50 = []
male_delivery_jux4_nhe50 = []
male_delivery_jux5_nhe50 = []
for seg in segment_early:
    file_sup = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    male_delivery_number_nhe50.append(0)
    male_delivery_sup_nhe50.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    male_delivery_jux1_nhe50.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    male_delivery_jux2_nhe50.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    male_delivery_jux3_nhe50.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    male_delivery_jux4_nhe50.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    male_delivery_jux5_nhe50.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        male_delivery_number_nhe50.append(0)
        male_delivery_sup_nhe50.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        male_delivery_jux1_nhe50.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        male_delivery_jux2_nhe50.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        male_delivery_jux3_nhe50.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        male_delivery_jux4_nhe50.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        male_delivery_jux5_nhe50.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(male_nhe50_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        male_delivery_number_nhe50.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    #===================================
    # NHE3 80% inhibited
    #===================================
female_delivery_number_nhe80 = []
female_delivery_sup_nhe80 = []
female_delivery_jux1_nhe80 = []
female_delivery_jux2_nhe80 = []
female_delivery_jux3_nhe80 = []
female_delivery_jux4_nhe80 = []
female_delivery_jux5_nhe80 = []
for seg in segment_early:
    file_sup = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    female_delivery_number_nhe80.append(0)
    female_delivery_sup_nhe80.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    female_delivery_jux1_nhe80.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    female_delivery_jux2_nhe80.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    female_delivery_jux3_nhe80.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    female_delivery_jux4_nhe80.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    female_delivery_jux5_nhe80.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        female_delivery_number_nhe80.append(0)
        female_delivery_sup_nhe80.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        female_delivery_jux1_nhe80.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        female_delivery_jux2_nhe80.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        female_delivery_jux3_nhe80.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        female_delivery_jux4_nhe80.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        female_delivery_jux5_nhe80.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(female_nhe80_file+'/female_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        female_delivery_number_nhe80.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]
    

male_delivery_number_nhe80 = []
male_delivery_sup_nhe80 = []
male_delivery_jux1_nhe80 = []
male_delivery_jux2_nhe80 = []
male_delivery_jux3_nhe80 = []
male_delivery_jux4_nhe80 = []
male_delivery_jux5_nhe80 = []
for seg in segment_early:
    file_sup = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_sup.txt','r')
    file_jux1 = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux1.txt','r')
    file_jux2 = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux2.txt','r')
    file_jux3 = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux3.txt','r')
    file_jux4 = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux4.txt','r')
    file_jux5 = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    male_delivery_number_nhe80.append(0)
    male_delivery_sup_nhe80.append(neph_weight[0]*datalist_sup[0]*solute_conversion)
    male_delivery_jux1_nhe80.append(neph_weight[1]*datalist_jux1[0]*solute_conversion)
    male_delivery_jux2_nhe80.append(neph_weight[2]*datalist_jux2[0]*solute_conversion)
    male_delivery_jux3_nhe80.append(neph_weight[3]*datalist_jux3[0]*solute_conversion)
    male_delivery_jux4_nhe80.append(neph_weight[4]*datalist_jux4[0]*solute_conversion)
    male_delivery_jux5_nhe80.append(neph_weight[5]*datalist_jux5[0]*solute_conversion)
    if seg == 'cnt':
        male_delivery_number_nhe80.append(0)
        male_delivery_sup_nhe80.append(neph_weight[0]*datalist_sup[-1]*solute_conversion)
        male_delivery_jux1_nhe80.append(neph_weight[1]*datalist_jux1[-1]*solute_conversion)
        male_delivery_jux2_nhe80.append(neph_weight[2]*datalist_jux2[-1]*solute_conversion)
        male_delivery_jux3_nhe80.append(neph_weight[3]*datalist_jux3[-1]*solute_conversion)
        male_delivery_jux4_nhe80.append(neph_weight[4]*datalist_jux4[-1]*solute_conversion)
        male_delivery_jux5_nhe80.append(neph_weight[5]*datalist_jux5[-1]*solute_conversion)
for seg in segment_late:
    file_data = open(male_nhe80_file+'/male_hum_'+seg+'_flow_of_'+s+'_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        male_delivery_number_nhe80.append(number_of_delivery*solute_conversion)
    else:
        number_of_delivery = datalist[0]


male_sup_base=axarr[2,0].bar(np.arange(len(segment[:6]))-bar_width,male_delivery_sup_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='blue',label='Men, Baseline')
male_jux_base=axarr[2,0].bar(np.arange(len(segment[:6]))-bar_width,[male_delivery_jux1_base[i]+male_delivery_jux2_base[i]+male_delivery_jux3_base[i]+male_delivery_jux4_base[i]+male_delivery_jux5_base[i] for i in range(len(male_delivery_sup_base))],bar_width,bottom=male_delivery_sup_base,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_base=axarr[2,0].bar(np.arange(len(segment))-bar_width,male_delivery_number_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='blue')

male_sup=axarr[2,0].bar(np.arange(len(segment[:6])),male_delivery_sup,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='royalblue',label='Men, ACEi')
male_jux=axarr[2,0].bar(np.arange(len(segment[:6])),[male_delivery_jux1[i]+male_delivery_jux2[i]+male_delivery_jux3[i]+male_delivery_jux4[i]+male_delivery_jux5[i] for i in range(len(male_delivery_sup))],bar_width,bottom=male_delivery_sup,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later=axarr[2,0].bar(np.arange(len(segment)),male_delivery_number,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='royalblue')

male_sup_nhe50=axarr[2,0].bar(np.arange(len(segment[:6]))+bar_width,male_delivery_sup_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='deepskyblue',label='Men, NHE3 PT alter')
male_jux_nhe50=axarr[2,0].bar(np.arange(len(segment[:6]))+bar_width,[male_delivery_jux1_nhe50[i]+male_delivery_jux2_nhe50[i]+male_delivery_jux3_nhe50[i]+male_delivery_jux4_nhe50[i]+male_delivery_jux5_nhe50[i] for i in range(len(male_delivery_sup_nhe50))],bar_width,bottom=male_delivery_sup_nhe50,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_nhe50=axarr[2,0].bar(np.arange(len(segment))+bar_width,male_delivery_number_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='deepskyblue')

male_sup_nhe80=axarr[2,0].bar(np.arange(len(segment[:6]))+2*bar_width,male_delivery_sup_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='paleturquoise',label='Men, distal alter')
male_jux_nhe80=axarr[2,0].bar(np.arange(len(segment[:6]))+2*bar_width,[male_delivery_jux1_nhe80[i]+male_delivery_jux2_nhe80[i]+male_delivery_jux3_nhe80[i]+male_delivery_jux4_nhe80[i]+male_delivery_jux5_nhe80[i] for i in range(len(male_delivery_sup_nhe80))],bar_width,bottom=male_delivery_sup_nhe80,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_nhe80=axarr[2,0].bar(np.arange(len(segment))+2*bar_width,male_delivery_number_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='paleturquoise')

female_sup_base=axarr[2,1].bar(np.arange(len(segment[:6]))-bar_width,female_delivery_sup_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='red',label='Women, Baseline')
female_jux_base=axarr[2,1].bar(np.arange(len(segment[:6]))-bar_width,[female_delivery_jux1_base[i]+female_delivery_jux2_base[i]+female_delivery_jux3_base[i]+female_delivery_jux4_base[i]+female_delivery_jux5_base[i] for i in range(len(female_delivery_sup_base))],bar_width,bottom=female_delivery_sup_base,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
female_later_base=axarr[2,1].bar(np.arange(len(segment))-bar_width,female_delivery_number_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='red')

Female_sup=axarr[2,1].bar(np.arange(len(segment[:6])),female_delivery_sup,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='magenta',label='Women, ACEi')
Female_jux=axarr[2,1].bar(np.arange(len(segment[:6])),[female_delivery_jux1[i]+female_delivery_jux2[i]+female_delivery_jux3[i]+female_delivery_jux4[i]+female_delivery_jux5[i] for i in range(len(female_delivery_sup))],bar_width,bottom=female_delivery_sup,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later=axarr[2,1].bar(np.arange(len(segment)),female_delivery_number,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='magenta')

Female_sup_nhe50=axarr[2,1].bar(np.arange(len(segment[:6]))+1*bar_width,female_delivery_sup_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='hotpink',label='Women, NHE3 PT alter')
Female_jux_nhe50=axarr[2,1].bar(np.arange(len(segment[:6]))+1*bar_width,[female_delivery_jux1_nhe50[i]+female_delivery_jux2_nhe50[i]+female_delivery_jux3_nhe50[i]+female_delivery_jux4_nhe50[i]+female_delivery_jux5_nhe50[i] for i in range(len(female_delivery_sup_nhe50))],bar_width,bottom=female_delivery_sup_nhe50,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later_nhe50=axarr[2,1].bar(np.arange(len(segment))+1*bar_width,female_delivery_number_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='hotpink')

Female_sup_nhe80=axarr[2,1].bar(np.arange(len(segment[:6]))+2*bar_width,female_delivery_sup_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='pink',label='Women, distal alter')
Female_jux_nhe80=axarr[2,1].bar(np.arange(len(segment[:6]))+2*bar_width,[female_delivery_jux1_nhe80[i]+female_delivery_jux2_nhe80[i]+female_delivery_jux3_nhe80[i]+female_delivery_jux4_nhe80[i]+female_delivery_jux5_nhe80[i] for i in range(len(female_delivery_sup_nhe80))],bar_width,bottom=female_delivery_sup_nhe80,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later_nhe80=axarr[2,1].bar(np.arange(len(segment))+2*bar_width,female_delivery_number_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='pink')

axarr[2,0].set_xticks(np.arange(len(segment))+1*bar_width)
axarr[2,0].set_xticklabels(segment,fontsize=40)
axarr[2,0].tick_params(axis='both',labelsize=40)
#axarr[2,0].set_ylim(0,150)
#ax.set_xlabel('Segment',fontsize=20)
axarr[2,0].set_ylabel('Cl$^-$ delivery (mol/Day)',fontsize=40)
#axarr[2,0].legend(fontsize=30,markerscale=30)
bar_width_ins = bar_width
#axins = inset_axes(ax,width=2,height=2,loc=7)
#Male_inset=axins.bar(np.arange(len(segment)),male_list,bar_width_ins,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='black',label='Male')
#Female_inset=axins.bar(np.arange(len(segment))+bar_width_ins,female_list,bar_width_ins,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white',label='Female')
#axins.set_xticks(np.arange(len(segment))+0.5*bar_width_ins)
#axins.set_xticklabels(segment,fontsize=20)
#axins.set_xlim(5-1.5*bar_width_ins,6+2*bar_width)
#axins.set_ylim(0,2)
#axins.tick_params(axis='both',labelsize=20)
#autolabel(Male,'left')
#autolabel(Female,'right')
axarr[2,1].set_xticks(np.arange(len(segment))+1*bar_width)
axarr[2,1].set_xticklabels(segment,fontsize=40)
axarr[2,1].tick_params(axis='both',labelsize=40)
#axarr[2,1].set_ylim(0,150)
#ax.set_xlabel('Segment',fontsize=20)
axarr[2,1].set_ylabel('Cl$^-$ delivery (mol/Day)',fontsize=40)
#axarr[2,1].legend(fontsize=30,markerscale=30)


bar_width_ins = bar_width
axins = inset_axes(axarr[2,0],width=3.5,height=3.5,loc=7)

male_sup_base_inset=axins.bar(np.arange(len(segment[:6]))-bar_width,male_delivery_sup_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='blue',label='Men, Baseline')
male_jux_base_inset=axins.bar(np.arange(len(segment[:6]))-bar_width,[male_delivery_jux1_base[i]+male_delivery_jux2_base[i]+male_delivery_jux3_base[i]+male_delivery_jux4_base[i]+male_delivery_jux5_base[i] for i in range(len(male_delivery_sup_base))],bar_width,bottom=male_delivery_sup_base,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_base_inset=axins.bar(np.arange(len(segment))-bar_width,male_delivery_number_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='blue')

male_sup_inset=axins.bar(np.arange(len(segment[:6])),male_delivery_sup,bar_width_ins,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='royalblue',label='Men, ACEi')
male_jux_inset=axins.bar(np.arange(len(segment[:6])),[male_delivery_jux1[i]+male_delivery_jux2[i]+male_delivery_jux3[i]+male_delivery_jux4[i]+male_delivery_jux5[i] for i in range(len(male_delivery_sup))],bar_width_ins,bottom=male_delivery_sup,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_inset=axins.bar(np.arange(len(segment)),male_delivery_number,bar_width_ins,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='royalblue')

male_sup_nhe50_inset=axins.bar(np.arange(len(segment[:6]))+bar_width,male_delivery_sup_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='deepskyblue',label='Men, ACEi PT')
male_jux_nhe50_inset=axins.bar(np.arange(len(segment[:6]))+bar_width,[male_delivery_jux1_nhe50[i]+male_delivery_jux2_nhe50[i]+male_delivery_jux3_nhe50[i]+male_delivery_jux4_nhe50[i]+male_delivery_jux5_nhe50[i] for i in range(len(male_delivery_sup_nhe50))],bar_width,bottom=male_delivery_sup_nhe50,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_nhe50_inset=axins.bar(np.arange(len(segment))+bar_width,male_delivery_number_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='deepskyblue')

male_sup_nhe80_inset=axins.bar(np.arange(len(segment[:6]))+2*bar_width,male_delivery_sup_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='paleturquoise',label='Men, ACEi distal')
male_jux_nhe80_inset=axins.bar(np.arange(len(segment[:6]))+2*bar_width,[male_delivery_jux1_nhe80[i]+male_delivery_jux2_nhe80[i]+male_delivery_jux3_nhe80[i]+male_delivery_jux4_nhe80[i]+male_delivery_jux5_nhe80[i] for i in range(len(male_delivery_sup_nhe80))],bar_width,bottom=male_delivery_sup_nhe80,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_nhe80_inset=axins.bar(np.arange(len(segment))+2*bar_width,male_delivery_number_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='paleturquoise')

axins.set_xticks(np.arange(len(segment))+1*bar_width_ins)
axins.set_xticklabels(segment,fontsize=40)
axins.set_xlim(5-2.5*bar_width_ins,6+3*bar_width_ins)
axins.set_ylim(0,1)
axins.tick_params(axis='both',labelsize=40)

bar_width_ins = bar_width
axins = inset_axes(axarr[2,1],width=3.5,height=3.5,loc=7)

female_sup_base_inset=axins.bar(np.arange(len(segment[:6]))-bar_width,female_delivery_sup_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='red',label='Women, Baseline')
female_jux_base_inset=axins.bar(np.arange(len(segment[:6]))-bar_width,[female_delivery_jux1_base[i]+female_delivery_jux2_base[i]+female_delivery_jux3_base[i]+female_delivery_jux4_base[i]+female_delivery_jux5_base[i] for i in range(len(female_delivery_sup_base))],bar_width,bottom=female_delivery_sup_base,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
female_later_base_inset=axins.bar(np.arange(len(segment))-bar_width,female_delivery_number_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='red')

Female_sup_inset=axins.bar(np.arange(len(segment[:6])),female_delivery_sup,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='magenta',label='Baseline')
Female_jux_inset=axins.bar(np.arange(len(segment[:6])),[female_delivery_jux1[i]+female_delivery_jux2[i]+female_delivery_jux3[i]+female_delivery_jux4[i]+female_delivery_jux5[i] for i in range(len(female_delivery_sup))],bar_width,bottom=female_delivery_sup,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later_inset=axins.bar(np.arange(len(segment)),female_delivery_number,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='magenta')

Female_sup_nhe50_inset=axins.bar(np.arange(len(segment[:6]))+1*bar_width,female_delivery_sup_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='hotpink',label='50% NHE3 inhib')
Female_jux_nhe50_inset=axins.bar(np.arange(len(segment[:6]))+1*bar_width,[female_delivery_jux1_nhe50[i]+female_delivery_jux2_nhe50[i]+female_delivery_jux3_nhe50[i]+female_delivery_jux4_nhe50[i]+female_delivery_jux5_nhe50[i] for i in range(len(female_delivery_sup_nhe50))],bar_width,bottom=female_delivery_sup_nhe50,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later_nhe50_inset=axins.bar(np.arange(len(segment))+1*bar_width,female_delivery_number_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='hotpink')

Female_sup_nhe80_inset=axins.bar(np.arange(len(segment[:6]))+2*bar_width,female_delivery_sup_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='pink',label='80% NHE3 inhib')
Female_jux_nhe80_inset=axins.bar(np.arange(len(segment[:6]))+2*bar_width,[female_delivery_jux1_nhe80[i]+female_delivery_jux2_nhe80[i]+female_delivery_jux3_nhe80[i]+female_delivery_jux4_nhe80[i]+female_delivery_jux5_nhe80[i] for i in range(len(female_delivery_sup_nhe80))],bar_width,bottom=female_delivery_sup_nhe80,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later_nhe80_inset=axins.bar(np.arange(len(segment))+2*bar_width,female_delivery_number_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='pink')

axins.set_xticks(np.arange(len(segment))+1*bar_width_ins)
axins.set_xticklabels(segment,fontsize=40)
axins.set_xlim(5-2.5*bar_width_ins,6+3*bar_width_ins)
axins.set_ylim(0,1)
axins.tick_params(axis='both',labelsize=40)



#=================================================================
# Water volume
#=================================================================

#==========================
# Baseline
#==========================

female_delivery_number_base = []
female_delivery_sup_base = []
female_delivery_jux1_base = []
female_delivery_jux2_base = []
female_delivery_jux3_base = []
female_delivery_jux4_base = []
female_delivery_jux5_base = []
for seg in segment_early:
    file_sup = open(female_base+'/female_hum_'+seg+'_water_volume_in_Lumen_sup.txt','r')
    file_jux1 = open(female_base+'/female_hum_'+seg+'_water_volume_in_Lumen_jux1.txt','r')
    file_jux2 = open(female_base+'/female_hum_'+seg+'_water_volume_in_Lumen_jux2.txt','r')
    file_jux3 = open(female_base+'/female_hum_'+seg+'_water_volume_in_Lumen_jux3.txt','r')
    file_jux4 = open(female_base+'/female_hum_'+seg+'_water_volume_in_Lumen_jux4.txt','r')
    file_jux5 = open(female_base+'/female_hum_'+seg+'_water_volume_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    female_delivery_number_base.append(0)
    female_delivery_sup_base.append(neph_weight[0]*datalist_sup[0]*volume_conversion)
    female_delivery_jux1_base.append(neph_weight[1]*datalist_jux1[0]*volume_conversion)
    female_delivery_jux2_base.append(neph_weight[2]*datalist_jux2[0]*volume_conversion)
    female_delivery_jux3_base.append(neph_weight[3]*datalist_jux3[0]*volume_conversion)
    female_delivery_jux4_base.append(neph_weight[4]*datalist_jux4[0]*volume_conversion)
    female_delivery_jux5_base.append(neph_weight[5]*datalist_jux5[0]*volume_conversion)
    if seg == 'cnt':
        female_delivery_number_base.append(0)
        female_delivery_sup_base.append(neph_weight[0]*datalist_sup[-1]*volume_conversion)
        female_delivery_jux1_base.append(neph_weight[1]*datalist_jux1[-1]*volume_conversion)
        female_delivery_jux2_base.append(neph_weight[2]*datalist_jux2[-1]*volume_conversion)
        female_delivery_jux3_base.append(neph_weight[3]*datalist_jux3[-1]*volume_conversion)
        female_delivery_jux4_base.append(neph_weight[4]*datalist_jux4[-1]*volume_conversion)
        female_delivery_jux5_base.append(neph_weight[5]*datalist_jux5[-1]*volume_conversion)
for seg in segment_late:
    file_data = open(female_base+'/female_hum_'+seg+'_water_volume_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        female_delivery_number_base.append(number_of_delivery*volume_conversion)
    else:
        number_of_delivery = datalist[0]
    

male_delivery_number_base = []
male_delivery_sup_base = []
male_delivery_jux1_base = []
male_delivery_jux2_base = []
male_delivery_jux3_base = []
male_delivery_jux4_base = []
male_delivery_jux5_base = []
for seg in segment_early:
    file_sup = open(male_base+'/male_hum_'+seg+'_water_volume_in_Lumen_sup.txt','r')
    file_jux1 = open(male_base+'/male_hum_'+seg+'_water_volume_in_Lumen_jux1.txt','r')
    file_jux2 = open(male_base+'/male_hum_'+seg+'_water_volume_in_Lumen_jux2.txt','r')
    file_jux3 = open(male_base+'/male_hum_'+seg+'_water_volume_in_Lumen_jux3.txt','r')
    file_jux4 = open(male_base+'/male_hum_'+seg+'_water_volume_in_Lumen_jux4.txt','r')
    file_jux5 = open(male_base+'/male_hum_'+seg+'_water_volume_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    male_delivery_number_base.append(0)
    male_delivery_sup_base.append(neph_weight[0]*datalist_sup[0]*volume_conversion)
    male_delivery_jux1_base.append(neph_weight[1]*datalist_jux1[0]*volume_conversion)
    male_delivery_jux2_base.append(neph_weight[2]*datalist_jux2[0]*volume_conversion)
    male_delivery_jux3_base.append(neph_weight[3]*datalist_jux3[0]*volume_conversion)
    male_delivery_jux4_base.append(neph_weight[4]*datalist_jux4[0]*volume_conversion)
    male_delivery_jux5_base.append(neph_weight[5]*datalist_jux5[0]*volume_conversion)
    if seg == 'cnt':
        male_delivery_number_base.append(0)
        male_delivery_sup_base.append(neph_weight[0]*datalist_sup[-1]*volume_conversion)
        male_delivery_jux1_base.append(neph_weight[1]*datalist_jux1[-1]*volume_conversion)
        male_delivery_jux2_base.append(neph_weight[2]*datalist_jux2[-1]*volume_conversion)
        male_delivery_jux3_base.append(neph_weight[3]*datalist_jux3[-1]*volume_conversion)
        male_delivery_jux4_base.append(neph_weight[4]*datalist_jux4[-1]*volume_conversion)
        male_delivery_jux5_base.append(neph_weight[5]*datalist_jux5[-1]*volume_conversion)
for seg in segment_late:
    file_data = open(male_base+'/male_hum_'+seg+'_water_volume_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        male_delivery_number_base.append(number_of_delivery*volume_conversion)
    else:
        number_of_delivery = datalist[0]

#==========================

female_delivery_number = []
female_delivery_sup = []
female_delivery_jux1 = []
female_delivery_jux2 = []
female_delivery_jux3 = []
female_delivery_jux4 = []
female_delivery_jux5 = []
for seg in segment_early:
    file_sup = open(female_normal_file+'/female_hum_'+seg+'_water_volume_in_Lumen_sup.txt','r')
    file_jux1 = open(female_normal_file+'/female_hum_'+seg+'_water_volume_in_Lumen_jux1.txt','r')
    file_jux2 = open(female_normal_file+'/female_hum_'+seg+'_water_volume_in_Lumen_jux2.txt','r')
    file_jux3 = open(female_normal_file+'/female_hum_'+seg+'_water_volume_in_Lumen_jux3.txt','r')
    file_jux4 = open(female_normal_file+'/female_hum_'+seg+'_water_volume_in_Lumen_jux4.txt','r')
    file_jux5 = open(female_normal_file+'/female_hum_'+seg+'_water_volume_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    female_delivery_number.append(0)
    female_delivery_sup.append(neph_weight[0]*datalist_sup[0]*volume_conversion)
    female_delivery_jux1.append(neph_weight[1]*datalist_jux1[0]*volume_conversion)
    female_delivery_jux2.append(neph_weight[2]*datalist_jux2[0]*volume_conversion)
    female_delivery_jux3.append(neph_weight[3]*datalist_jux3[0]*volume_conversion)
    female_delivery_jux4.append(neph_weight[4]*datalist_jux4[0]*volume_conversion)
    female_delivery_jux5.append(neph_weight[5]*datalist_jux5[0]*volume_conversion)
    if seg == 'cnt':
        female_delivery_number.append(0)
        female_delivery_sup.append(neph_weight[0]*datalist_sup[-1]*volume_conversion)
        female_delivery_jux1.append(neph_weight[1]*datalist_jux1[-1]*volume_conversion)
        female_delivery_jux2.append(neph_weight[2]*datalist_jux2[-1]*volume_conversion)
        female_delivery_jux3.append(neph_weight[3]*datalist_jux3[-1]*volume_conversion)
        female_delivery_jux4.append(neph_weight[4]*datalist_jux4[-1]*volume_conversion)
        female_delivery_jux5.append(neph_weight[5]*datalist_jux5[-1]*volume_conversion)
for seg in segment_late:
    file_data = open(female_normal_file+'/female_hum_'+seg+'_water_volume_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        female_delivery_number.append(number_of_delivery*volume_conversion)
    else:
        number_of_delivery = datalist[0]
    

male_delivery_number = []
male_delivery_sup = []
male_delivery_jux1 = []
male_delivery_jux2 = []
male_delivery_jux3 = []
male_delivery_jux4 = []
male_delivery_jux5 = []
for seg in segment_early:
    file_sup = open(male_normal_file+'/male_hum_'+seg+'_water_volume_in_Lumen_sup.txt','r')
    file_jux1 = open(male_normal_file+'/male_hum_'+seg+'_water_volume_in_Lumen_jux1.txt','r')
    file_jux2 = open(male_normal_file+'/male_hum_'+seg+'_water_volume_in_Lumen_jux2.txt','r')
    file_jux3 = open(male_normal_file+'/male_hum_'+seg+'_water_volume_in_Lumen_jux3.txt','r')
    file_jux4 = open(male_normal_file+'/male_hum_'+seg+'_water_volume_in_Lumen_jux4.txt','r')
    file_jux5 = open(male_normal_file+'/male_hum_'+seg+'_water_volume_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    male_delivery_number.append(0)
    male_delivery_sup.append(neph_weight[0]*datalist_sup[0]*volume_conversion)
    male_delivery_jux1.append(neph_weight[1]*datalist_jux1[0]*volume_conversion)
    male_delivery_jux2.append(neph_weight[2]*datalist_jux2[0]*volume_conversion)
    male_delivery_jux3.append(neph_weight[3]*datalist_jux3[0]*volume_conversion)
    male_delivery_jux4.append(neph_weight[4]*datalist_jux4[0]*volume_conversion)
    male_delivery_jux5.append(neph_weight[5]*datalist_jux5[0]*volume_conversion)
    if seg == 'cnt':
        male_delivery_number.append(0)
        male_delivery_sup.append(neph_weight[0]*datalist_sup[-1]*volume_conversion)
        male_delivery_jux1.append(neph_weight[1]*datalist_jux1[-1]*volume_conversion)
        male_delivery_jux2.append(neph_weight[2]*datalist_jux2[-1]*volume_conversion)
        male_delivery_jux3.append(neph_weight[3]*datalist_jux3[-1]*volume_conversion)
        male_delivery_jux4.append(neph_weight[4]*datalist_jux4[-1]*volume_conversion)
        male_delivery_jux5.append(neph_weight[5]*datalist_jux5[-1]*volume_conversion)
for seg in segment_late:
    file_data = open(male_normal_file+'/male_hum_'+seg+'_water_volume_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        male_delivery_number.append(number_of_delivery*volume_conversion)
    else:
        number_of_delivery = datalist[0]
    
    #===================================
    # NHE3 50% inhibited
    #===================================
female_delivery_number_nhe50 = []
female_delivery_sup_nhe50 = []
female_delivery_jux1_nhe50 = []
female_delivery_jux2_nhe50 = []
female_delivery_jux3_nhe50 = []
female_delivery_jux4_nhe50 = []
female_delivery_jux5_nhe50 = []
for seg in segment_early:
    file_sup = open(female_nhe50_file+'/female_hum_'+seg+'_water_volume_in_Lumen_sup.txt','r')
    file_jux1 = open(female_nhe50_file+'/female_hum_'+seg+'_water_volume_in_Lumen_jux1.txt','r')
    file_jux2 = open(female_nhe50_file+'/female_hum_'+seg+'_water_volume_in_Lumen_jux2.txt','r')
    file_jux3 = open(female_nhe50_file+'/female_hum_'+seg+'_water_volume_in_Lumen_jux3.txt','r')
    file_jux4 = open(female_nhe50_file+'/female_hum_'+seg+'_water_volume_in_Lumen_jux4.txt','r')
    file_jux5 = open(female_nhe50_file+'/female_hum_'+seg+'_water_volume_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    female_delivery_number_nhe50.append(0)
    female_delivery_sup_nhe50.append(neph_weight[0]*datalist_sup[0]*volume_conversion)
    female_delivery_jux1_nhe50.append(neph_weight[1]*datalist_jux1[0]*volume_conversion)
    female_delivery_jux2_nhe50.append(neph_weight[2]*datalist_jux2[0]*volume_conversion)
    female_delivery_jux3_nhe50.append(neph_weight[3]*datalist_jux3[0]*volume_conversion)
    female_delivery_jux4_nhe50.append(neph_weight[4]*datalist_jux4[0]*volume_conversion)
    female_delivery_jux5_nhe50.append(neph_weight[5]*datalist_jux5[0]*volume_conversion)
    if seg == 'cnt':
        female_delivery_number_nhe50.append(0)
        female_delivery_sup_nhe50.append(neph_weight[0]*datalist_sup[-1]*volume_conversion)
        female_delivery_jux1_nhe50.append(neph_weight[1]*datalist_jux1[-1]*volume_conversion)
        female_delivery_jux2_nhe50.append(neph_weight[2]*datalist_jux2[-1]*volume_conversion)
        female_delivery_jux3_nhe50.append(neph_weight[3]*datalist_jux3[-1]*volume_conversion)
        female_delivery_jux4_nhe50.append(neph_weight[4]*datalist_jux4[-1]*volume_conversion)
        female_delivery_jux5_nhe50.append(neph_weight[5]*datalist_jux5[-1]*volume_conversion)
for seg in segment_late:
    file_data = open(female_nhe50_file+'/female_hum_'+seg+'_water_volume_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        female_delivery_number_nhe50.append(number_of_delivery*volume_conversion)
    else:
        number_of_delivery = datalist[0]
    

male_delivery_number_nhe50 = []
male_delivery_sup_nhe50 = []
male_delivery_jux1_nhe50 = []
male_delivery_jux2_nhe50 = []
male_delivery_jux3_nhe50 = []
male_delivery_jux4_nhe50 = []
male_delivery_jux5_nhe50 = []
for seg in segment_early:
    file_sup = open(male_nhe50_file+'/male_hum_'+seg+'_water_volume_in_Lumen_sup.txt','r')
    file_jux1 = open(male_nhe50_file+'/male_hum_'+seg+'_water_volume_in_Lumen_jux1.txt','r')
    file_jux2 = open(male_nhe50_file+'/male_hum_'+seg+'_water_volume_in_Lumen_jux2.txt','r')
    file_jux3 = open(male_nhe50_file+'/male_hum_'+seg+'_water_volume_in_Lumen_jux3.txt','r')
    file_jux4 = open(male_nhe50_file+'/male_hum_'+seg+'_water_volume_in_Lumen_jux4.txt','r')
    file_jux5 = open(male_nhe50_file+'/male_hum_'+seg+'_water_volume_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    male_delivery_number_nhe50.append(0)
    male_delivery_sup_nhe50.append(neph_weight[0]*datalist_sup[0]*volume_conversion)
    male_delivery_jux1_nhe50.append(neph_weight[1]*datalist_jux1[0]*volume_conversion)
    male_delivery_jux2_nhe50.append(neph_weight[2]*datalist_jux2[0]*volume_conversion)
    male_delivery_jux3_nhe50.append(neph_weight[3]*datalist_jux3[0]*volume_conversion)
    male_delivery_jux4_nhe50.append(neph_weight[4]*datalist_jux4[0]*volume_conversion)
    male_delivery_jux5_nhe50.append(neph_weight[5]*datalist_jux5[0]*volume_conversion)
    if seg == 'cnt':
        male_delivery_number_nhe50.append(0)
        male_delivery_sup_nhe50.append(neph_weight[0]*datalist_sup[-1]*volume_conversion)
        male_delivery_jux1_nhe50.append(neph_weight[1]*datalist_jux1[-1]*volume_conversion)
        male_delivery_jux2_nhe50.append(neph_weight[2]*datalist_jux2[-1]*volume_conversion)
        male_delivery_jux3_nhe50.append(neph_weight[3]*datalist_jux3[-1]*volume_conversion)
        male_delivery_jux4_nhe50.append(neph_weight[4]*datalist_jux4[-1]*volume_conversion)
        male_delivery_jux5_nhe50.append(neph_weight[5]*datalist_jux5[-1]*volume_conversion)
for seg in segment_late:
    file_data = open(male_nhe50_file+'/male_hum_'+seg+'_water_volume_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        male_delivery_number_nhe50.append(number_of_delivery*volume_conversion)
    else:
        number_of_delivery = datalist[0]
    #===================================
    # NHE3 80% inhibited
    #===================================
female_delivery_number_nhe80 = []
female_delivery_sup_nhe80 = []
female_delivery_jux1_nhe80 = []
female_delivery_jux2_nhe80 = []
female_delivery_jux3_nhe80 = []
female_delivery_jux4_nhe80 = []
female_delivery_jux5_nhe80 = []
for seg in segment_early:
    file_sup = open(female_nhe80_file+'/female_hum_'+seg+'_water_volume_in_Lumen_sup.txt','r')
    file_jux1 = open(female_nhe80_file+'/female_hum_'+seg+'_water_volume_in_Lumen_jux1.txt','r')
    file_jux2 = open(female_nhe80_file+'/female_hum_'+seg+'_water_volume_in_Lumen_jux2.txt','r')
    file_jux3 = open(female_nhe80_file+'/female_hum_'+seg+'_water_volume_in_Lumen_jux3.txt','r')
    file_jux4 = open(female_nhe80_file+'/female_hum_'+seg+'_water_volume_in_Lumen_jux4.txt','r')
    file_jux5 = open(female_nhe80_file+'/female_hum_'+seg+'_water_volume_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    female_delivery_number_nhe80.append(0)
    female_delivery_sup_nhe80.append(neph_weight[0]*datalist_sup[0]*volume_conversion)
    female_delivery_jux1_nhe80.append(neph_weight[1]*datalist_jux1[0]*volume_conversion)
    female_delivery_jux2_nhe80.append(neph_weight[2]*datalist_jux2[0]*volume_conversion)
    female_delivery_jux3_nhe80.append(neph_weight[3]*datalist_jux3[0]*volume_conversion)
    female_delivery_jux4_nhe80.append(neph_weight[4]*datalist_jux4[0]*volume_conversion)
    female_delivery_jux5_nhe80.append(neph_weight[5]*datalist_jux5[0]*volume_conversion)
    if seg == 'cnt':
        female_delivery_number_nhe80.append(0)
        female_delivery_sup_nhe80.append(neph_weight[0]*datalist_sup[-1]*volume_conversion)
        female_delivery_jux1_nhe80.append(neph_weight[1]*datalist_jux1[-1]*volume_conversion)
        female_delivery_jux2_nhe80.append(neph_weight[2]*datalist_jux2[-1]*volume_conversion)
        female_delivery_jux3_nhe80.append(neph_weight[3]*datalist_jux3[-1]*volume_conversion)
        female_delivery_jux4_nhe80.append(neph_weight[4]*datalist_jux4[-1]*volume_conversion)
        female_delivery_jux5_nhe80.append(neph_weight[5]*datalist_jux5[-1]*volume_conversion)
for seg in segment_late:
    file_data = open(female_nhe80_file+'/female_hum_'+seg+'_water_volume_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        female_delivery_number_nhe80.append(number_of_delivery*volume_conversion)
    else:
        number_of_delivery = datalist[0]
    

male_delivery_number_nhe80 = []
male_delivery_sup_nhe80 = []
male_delivery_jux1_nhe80 = []
male_delivery_jux2_nhe80 = []
male_delivery_jux3_nhe80 = []
male_delivery_jux4_nhe80 = []
male_delivery_jux5_nhe80 = []
for seg in segment_early:
    file_sup = open(male_nhe80_file+'/male_hum_'+seg+'_water_volume_in_Lumen_sup.txt','r')
    file_jux1 = open(male_nhe80_file+'/male_hum_'+seg+'_water_volume_in_Lumen_jux1.txt','r')
    file_jux2 = open(male_nhe80_file+'/male_hum_'+seg+'_water_volume_in_Lumen_jux2.txt','r')
    file_jux3 = open(male_nhe80_file+'/male_hum_'+seg+'_water_volume_in_Lumen_jux3.txt','r')
    file_jux4 = open(male_nhe80_file+'/male_hum_'+seg+'_water_volume_in_Lumen_jux4.txt','r')
    file_jux5 = open(male_nhe80_file+'/male_hum_'+seg+'_water_volume_in_Lumen_jux5.txt','r')
    datalist_sup = []
    datalist_jux1 = []
    datalist_jux2 = []
    datalist_jux3 = []
    datalist_jux4 = []
    datalist_jux5 = []
    for i in file_sup:
        line = i.split(' ')
        datalist_sup.append(float(line[0]))
    for i in file_jux1:
        line = i.split(' ')
        datalist_jux1.append(float(line[0]))
    for i in file_jux2:
        line = i.split(' ')
        datalist_jux2.append(float(line[0]))
    for i in file_jux3:
        line = i.split(' ')
        datalist_jux3.append(float(line[0]))
    for i in file_jux4:
        line = i.split(' ')
        datalist_jux4.append(float(line[0]))
    for i in file_jux5:
        line = i.split(' ')
        datalist_jux5.append(float(line[0]))
    number_of_delivery = neph_weight[0]*datalist_sup[0]+neph_weight[1]*datalist_jux1[0]+neph_weight[2]*datalist_jux2[0]+neph_weight[3]*datalist_jux3[0]+neph_weight[4]*datalist_jux4[0]+neph_weight[5]*datalist_jux5[0]
    male_delivery_number_nhe80.append(0)
    male_delivery_sup_nhe80.append(neph_weight[0]*datalist_sup[0]*volume_conversion)
    male_delivery_jux1_nhe80.append(neph_weight[1]*datalist_jux1[0]*volume_conversion)
    male_delivery_jux2_nhe80.append(neph_weight[2]*datalist_jux2[0]*volume_conversion)
    male_delivery_jux3_nhe80.append(neph_weight[3]*datalist_jux3[0]*volume_conversion)
    male_delivery_jux4_nhe80.append(neph_weight[4]*datalist_jux4[0]*volume_conversion)
    male_delivery_jux5_nhe80.append(neph_weight[5]*datalist_jux5[0]*volume_conversion)
    if seg == 'cnt':
        male_delivery_number_nhe80.append(0)
        male_delivery_sup_nhe80.append(neph_weight[0]*datalist_sup[-1]*volume_conversion)
        male_delivery_jux1_nhe80.append(neph_weight[1]*datalist_jux1[-1]*volume_conversion)
        male_delivery_jux2_nhe80.append(neph_weight[2]*datalist_jux2[-1]*volume_conversion)
        male_delivery_jux3_nhe80.append(neph_weight[3]*datalist_jux3[-1]*volume_conversion)
        male_delivery_jux4_nhe80.append(neph_weight[4]*datalist_jux4[-1]*volume_conversion)
        male_delivery_jux5_nhe80.append(neph_weight[5]*datalist_jux5[-1]*volume_conversion)
for seg in segment_late:
    file_data = open(male_nhe80_file+'/male_hum_'+seg+'_water_volume_in_Lumen.txt','r')
    datalist = []
    for i in file_data:
        line = i.split(' ')
        datalist.append(float(line[0]))
    if seg == 'imcd':
        number_of_delivery = datalist[-1]
        male_delivery_number_nhe80.append(number_of_delivery*volume_conversion)
    else:
        number_of_delivery = datalist[0]

male_sup_base=axarr[3,0].bar(np.arange(len(segment[:6]))-bar_width,male_delivery_sup_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='blue',label='Men, Baseline')
male_jux_base=axarr[3,0].bar(np.arange(len(segment[:6]))-bar_width,[male_delivery_jux1_base[i]+male_delivery_jux2_base[i]+male_delivery_jux3_base[i]+male_delivery_jux4_base[i]+male_delivery_jux5_base[i] for i in range(len(male_delivery_sup_base))],bar_width,bottom=male_delivery_sup_base,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_base=axarr[3,0].bar(np.arange(len(segment))-bar_width,male_delivery_number_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='blue')

male_sup=axarr[3,0].bar(np.arange(len(segment[:6])),male_delivery_sup,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='royalblue',label='Men, ACEi')
male_jux=axarr[3,0].bar(np.arange(len(segment[:6])),[male_delivery_jux1[i]+male_delivery_jux2[i]+male_delivery_jux3[i]+male_delivery_jux4[i]+male_delivery_jux5[i] for i in range(len(male_delivery_sup))],bar_width,bottom=male_delivery_sup,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later=axarr[3,0].bar(np.arange(len(segment)),male_delivery_number,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='royalblue')

male_sup_nhe50=axarr[3,0].bar(np.arange(len(segment[:6]))+bar_width,male_delivery_sup_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='deepskyblue',label='Men, NHE3 PT alter')
male_jux_nhe50=axarr[3,0].bar(np.arange(len(segment[:6]))+bar_width,[male_delivery_jux1_nhe50[i]+male_delivery_jux2_nhe50[i]+male_delivery_jux3_nhe50[i]+male_delivery_jux4_nhe50[i]+male_delivery_jux5_nhe50[i] for i in range(len(male_delivery_sup_nhe50))],bar_width,bottom=male_delivery_sup_nhe50,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_nhe50=axarr[3,0].bar(np.arange(len(segment))+bar_width,male_delivery_number_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='deepskyblue')

male_sup_nhe80=axarr[3,0].bar(np.arange(len(segment[:6]))+2*bar_width,male_delivery_sup_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='paleturquoise',label='Men, distal alter')
male_jux_nhe80=axarr[3,0].bar(np.arange(len(segment[:6]))+2*bar_width,[male_delivery_jux1_nhe80[i]+male_delivery_jux2_nhe80[i]+male_delivery_jux3_nhe80[i]+male_delivery_jux4_nhe80[i]+male_delivery_jux5_nhe80[i] for i in range(len(male_delivery_sup_nhe80))],bar_width,bottom=male_delivery_sup_nhe80,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_nhe80=axarr[3,0].bar(np.arange(len(segment))+2*bar_width,male_delivery_number_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='paleturquoise')

female_sup_base=axarr[3,1].bar(np.arange(len(segment[:6]))-bar_width,female_delivery_sup_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='red',label='Women, Baseline')
female_jux_base=axarr[3,1].bar(np.arange(len(segment[:6]))-bar_width,[female_delivery_jux1_base[i]+female_delivery_jux2_base[i]+female_delivery_jux3_base[i]+female_delivery_jux4_base[i]+female_delivery_jux5_base[i] for i in range(len(female_delivery_sup_base))],bar_width,bottom=female_delivery_sup_base,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
female_later_base=axarr[3,1].bar(np.arange(len(segment))-bar_width,female_delivery_number_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='red')

Female_sup=axarr[3,1].bar(np.arange(len(segment[:6])),female_delivery_sup,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='magenta',label='Women, ACEi')
Female_jux=axarr[3,1].bar(np.arange(len(segment[:6])),[female_delivery_jux1[i]+female_delivery_jux2[i]+female_delivery_jux3[i]+female_delivery_jux4[i]+female_delivery_jux5[i] for i in range(len(female_delivery_sup))],bar_width,bottom=female_delivery_sup,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later=axarr[3,1].bar(np.arange(len(segment)),female_delivery_number,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='magenta')

Female_sup_nhe50=axarr[3,1].bar(np.arange(len(segment[:6]))+1*bar_width,female_delivery_sup_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='hotpink',label='Women, NHE3 PT alter')
Female_jux_nhe50=axarr[3,1].bar(np.arange(len(segment[:6]))+1*bar_width,[female_delivery_jux1_nhe50[i]+female_delivery_jux2_nhe50[i]+female_delivery_jux3_nhe50[i]+female_delivery_jux4_nhe50[i]+female_delivery_jux5_nhe50[i] for i in range(len(female_delivery_sup_nhe50))],bar_width,bottom=female_delivery_sup_nhe50,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later_nhe50=axarr[3,1].bar(np.arange(len(segment))+1*bar_width,female_delivery_number_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='hotpink')

Female_sup_nhe80=axarr[3,1].bar(np.arange(len(segment[:6]))+2*bar_width,female_delivery_sup_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='pink',label='Women, distal alter')
Female_jux_nhe80=axarr[3,1].bar(np.arange(len(segment[:6]))+2*bar_width,[female_delivery_jux1_nhe80[i]+female_delivery_jux2_nhe80[i]+female_delivery_jux3_nhe80[i]+female_delivery_jux4_nhe80[i]+female_delivery_jux5_nhe80[i] for i in range(len(female_delivery_sup_nhe80))],bar_width,bottom=female_delivery_sup_nhe80,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later_nhe80=axarr[3,1].bar(np.arange(len(segment))+2*bar_width,female_delivery_number_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='pink')

axarr[3,0].set_xticks(np.arange(len(segment))+1*bar_width)
axarr[3,0].set_xticklabels(segment,fontsize=40)
axarr[3,0].tick_params(axis='both',labelsize=40)
#axarr[3,0].set_ylim(0,1.5)
#ax.set_xlabel('Segment',fontsize=20)
axarr[3,0].set_ylabel('Volume delivery (mol/Day)',fontsize=40)
#axarr[3,0].legend(fontsize=30,markerscale=30)
bar_width_ins = bar_width
#axins = inset_axes(ax,width=2,height=2,loc=7)
#Male_inset=axins.bar(np.arange(len(segment)),male_list,bar_width_ins,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='black',label='Male')
#Female_inset=axins.bar(np.arange(len(segment))+bar_width_ins,female_list,bar_width_ins,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white',label='Female')
#axins.set_xticks(np.arange(len(segment))+0.5*bar_width_ins)
#axins.set_xticklabels(segment,fontsize=20)
#axins.set_xlim(5-1.5*bar_width_ins,6+2*bar_width)
#axins.set_ylim(0,2)
#axins.tick_params(axis='both',labelsize=20)
#autolabel(Male,'left')
#autolabel(Female,'right')
axarr[3,1].set_xticks(np.arange(len(segment))+1*bar_width)
axarr[3,1].set_xticklabels(segment,fontsize=40)
axarr[3,1].tick_params(axis='both',labelsize=40)
#axarr[3,1].set_ylim(0,1.5)
#ax.set_xlabel('Segment',fontsize=20)
axarr[3,1].set_ylabel('Volume delivery (L/Day)',fontsize=40)
#axarr[3,1].legend(fontsize=30,markerscale=30)

bar_width_ins = bar_width
axins = inset_axes(axarr[3,0],width=3.5,height=3.5,loc=7)

male_sup_base_inset=axins.bar(np.arange(len(segment[:6]))-bar_width,male_delivery_sup_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='blue',label='Men, Baseline')
male_jux_base_inset=axins.bar(np.arange(len(segment[:6]))-bar_width,[male_delivery_jux1_base[i]+male_delivery_jux2_base[i]+male_delivery_jux3_base[i]+male_delivery_jux4_base[i]+male_delivery_jux5_base[i] for i in range(len(male_delivery_sup_base))],bar_width,bottom=male_delivery_sup_base,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_base_inset=axins.bar(np.arange(len(segment))-bar_width,male_delivery_number_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='blue')

male_sup_inset=axins.bar(np.arange(len(segment[:6])),male_delivery_sup,bar_width_ins,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='royalblue',label='Men, ACEi')
male_jux_inset=axins.bar(np.arange(len(segment[:6])),[male_delivery_jux1[i]+male_delivery_jux2[i]+male_delivery_jux3[i]+male_delivery_jux4[i]+male_delivery_jux5[i] for i in range(len(male_delivery_sup))],bar_width_ins,bottom=male_delivery_sup,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_inset=axins.bar(np.arange(len(segment)),male_delivery_number,bar_width_ins,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='royalblue')

male_sup_nhe50_inset=axins.bar(np.arange(len(segment[:6]))+bar_width,male_delivery_sup_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='deepskyblue',label='Men, ACEi PT')
male_jux_nhe50_inset=axins.bar(np.arange(len(segment[:6]))+bar_width,[male_delivery_jux1_nhe50[i]+male_delivery_jux2_nhe50[i]+male_delivery_jux3_nhe50[i]+male_delivery_jux4_nhe50[i]+male_delivery_jux5_nhe50[i] for i in range(len(male_delivery_sup_nhe50))],bar_width,bottom=male_delivery_sup_nhe50,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_nhe50_inset=axins.bar(np.arange(len(segment))+bar_width,male_delivery_number_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='deepskyblue')

male_sup_nhe80_inset=axins.bar(np.arange(len(segment[:6]))+2*bar_width,male_delivery_sup_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='paleturquoise',label='Men, ACEi distal')
male_jux_nhe80_inset=axins.bar(np.arange(len(segment[:6]))+2*bar_width,[male_delivery_jux1_nhe80[i]+male_delivery_jux2_nhe80[i]+male_delivery_jux3_nhe80[i]+male_delivery_jux4_nhe80[i]+male_delivery_jux5_nhe80[i] for i in range(len(male_delivery_sup_nhe80))],bar_width,bottom=male_delivery_sup_nhe80,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
male_later_nhe80_inset=axins.bar(np.arange(len(segment))+2*bar_width,male_delivery_number_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='paleturquoise')

axins.set_xticks(np.arange(len(segment))+1*bar_width_ins)
axins.set_xticklabels(segment,fontsize=40)
axins.set_xlim(5-2.5*bar_width_ins,6+3*bar_width_ins)
axins.set_ylim(0,15)
axins.tick_params(axis='both',labelsize=40)

bar_width_ins = bar_width
axins = inset_axes(axarr[3,1],width=3.5,height=3.5,loc=7)

female_sup_base_inset=axins.bar(np.arange(len(segment[:6]))-bar_width,female_delivery_sup_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='red',label='Women, Baseline')
female_jux_base_inset=axins.bar(np.arange(len(segment[:6]))-bar_width,[female_delivery_jux1_base[i]+female_delivery_jux2_base[i]+female_delivery_jux3_base[i]+female_delivery_jux4_base[i]+female_delivery_jux5_base[i] for i in range(len(female_delivery_sup_base))],bar_width,bottom=female_delivery_sup_base,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
female_later_base_inset=axins.bar(np.arange(len(segment))-bar_width,female_delivery_number_base,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='red')

Female_sup_inset=axins.bar(np.arange(len(segment[:6])),female_delivery_sup,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='magenta',label='Baseline')
Female_jux_inset=axins.bar(np.arange(len(segment[:6])),[female_delivery_jux1[i]+female_delivery_jux2[i]+female_delivery_jux3[i]+female_delivery_jux4[i]+female_delivery_jux5[i] for i in range(len(female_delivery_sup))],bar_width,bottom=female_delivery_sup,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later_inset=axins.bar(np.arange(len(segment)),female_delivery_number,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='magenta')

Female_sup_nhe50_inset=axins.bar(np.arange(len(segment[:6]))+1*bar_width,female_delivery_sup_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='hotpink',label='50% NHE3 inhib')
Female_jux_nhe50_inset=axins.bar(np.arange(len(segment[:6]))+1*bar_width,[female_delivery_jux1_nhe50[i]+female_delivery_jux2_nhe50[i]+female_delivery_jux3_nhe50[i]+female_delivery_jux4_nhe50[i]+female_delivery_jux5_nhe50[i] for i in range(len(female_delivery_sup_nhe50))],bar_width,bottom=female_delivery_sup_nhe50,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later_nhe50_inset=axins.bar(np.arange(len(segment))+1*bar_width,female_delivery_number_nhe50,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='hotpink')

Female_sup_nhe80_inset=axins.bar(np.arange(len(segment[:6]))+2*bar_width,female_delivery_sup_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='pink',label='80% NHE3 inhib')
Female_jux_nhe80_inset=axins.bar(np.arange(len(segment[:6]))+2*bar_width,[female_delivery_jux1_nhe80[i]+female_delivery_jux2_nhe80[i]+female_delivery_jux3_nhe80[i]+female_delivery_jux4_nhe80[i]+female_delivery_jux5_nhe80[i] for i in range(len(female_delivery_sup_nhe80))],bar_width,bottom=female_delivery_sup_nhe80,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='white')
Female_later_nhe80_inset=axins.bar(np.arange(len(segment))+2*bar_width,female_delivery_number_nhe80,bar_width,align='center',alpha=0.8,linewidth=2,edgecolor='black',color='pink')

axins.set_xticks(np.arange(len(segment))+1*bar_width_ins)
axins.set_xticklabels(segment,fontsize=40)
axins.set_xlim(5-2.5*bar_width_ins,6+3*bar_width_ins)
axins.set_ylim(0,15)
axins.tick_params(axis='both',labelsize=40)

axarr[0,0].text(-1.5,axarr[0,0].get_ylim()[1],'A',size=40,weight='bold')
axarr[0,1].text(-1.5,axarr[0,1].get_ylim()[1],'B',size=40,weight='bold')
axarr[1,0].text(-1.5,axarr[1,0].get_ylim()[1],'C',size=40,weight='bold')
axarr[1,1].text(-1.5,axarr[1,1].get_ylim()[1],'D',size=40,weight='bold')
axarr[2,0].text(-1.5,axarr[2,0].get_ylim()[1],'E',size=40,weight='bold')
axarr[2,1].text(-1.5,axarr[2,1].get_ylim()[1],'F',size=40,weight='bold')
axarr[3,0].text(-1.5,axarr[3,0].get_ylim()[1],'G',size=40,weight='bold')
axarr[3,1].text(-1.5,axarr[3,1].get_ylim()[1],'H',size=40,weight='bold')

plt.savefig('ACEi',bbox_inches='tight')