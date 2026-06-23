import matplotlib.pyplot as plt
import skrf as rf
import numpy as np

# modal conversion of the S-parameter
def s2smm(sedata, flag_port1to3_2to4_thru):

    if flag_port1to3_2to4_thru == True:
        # If the port order is as presented below (PORT1 and PORT2 are on the same side, PORT3 and PORT4 are on the other side), re-order your ports (set flag_port1to3_2to4_thru = True) 
        # TX -> RX
        # PORT1 -> PORT3
        # PORT2 -> PORT4
        # network_se.renumber([0, 1, 2, 3], [new_order_map])
        sedata.renumber([0, 1, 2, 3], [0, 2, 1, 3])  # convert 1->3 and 2->4, to 1->2 and 3->4 

    # before the conversion, make sure that the S-parameters's port order is correct as shown below
    # TX -> RX
    # PORT1 -> PORT2
    # PORT3 -> PORT4
    sdd = np.zeros((len(sedata.frequency), 2, 2), dtype=complex)
    sdd[:,0,0]= 0.5*(sedata.s[:,0,0] - sedata.s[:,0,2] - sedata.s[:,2,0] + sedata.s[:,2,2]) #sdd11
    sdd[:,0,1]= 0.5*(sedata.s[:,0,1] - sedata.s[:,0,3] - sedata.s[:,2,1] + sedata.s[:,2,3]) #sdd12
    sdd[:,1,0]= 0.5*(sedata.s[:,1,0] - sedata.s[:,1,2] - sedata.s[:,3,0] + sedata.s[:,3,2]) #sdd21
    sdd[:,1,1]= 0.5*(sedata.s[:,1,1] - sedata.s[:,1,3] - sedata.s[:,3,1] + sedata.s[:,3,3]) #sdd22
    #sdd_ntw = rf.Network(frequency=freq, s=sdd, name='SDD')

    scc = np.zeros((len(sedata.frequency), 2, 2), dtype=complex)
    scc[:,0,0]= 0.5*(sedata.s[:,0,0] + sedata.s[:,0,2] + sedata.s[:,2,0] + sedata.s[:,2,2]) #scc11
    scc[:,0,1]= 0.5*(sedata.s[:,0,1] + sedata.s[:,0,3] + sedata.s[:,2,1] + sedata.s[:,2,3]) #scc12
    scc[:,1,0]= 0.5*(sedata.s[:,1,0] + sedata.s[:,1,2] + sedata.s[:,3,0] + sedata.s[:,3,2]) #scc21
    scc[:,1,1]= 0.5*(sedata.s[:,1,1] + sedata.s[:,1,3] + sedata.s[:,3,1] + sedata.s[:,3,3]) #scc22

    sdc = np.zeros((len(sedata.frequency), 2, 2), dtype=complex)
    sdc[:,0,0]= 0.5*(sedata.s[:,0,0] + sedata.s[:,0,2] - sedata.s[:,2,0] - sedata.s[:,2,2]) #sdc11
    sdc[:,0,1]= 0.5*(sedata.s[:,0,1] + sedata.s[:,0,3] - sedata.s[:,2,1] - sedata.s[:,2,3]) #sdc12
    sdc[:,1,0]= 0.5*(sedata.s[:,1,0] + sedata.s[:,1,2] - sedata.s[:,3,0] - sedata.s[:,3,2]) #sdc21
    sdc[:,1,1]= 0.5*(sedata.s[:,1,1] + sedata.s[:,1,3] - sedata.s[:,3,1] - sedata.s[:,3,3]) #sdc22

    scd = np.zeros((len(sedata.frequency), 2, 2), dtype=complex)
    scd[:,0,0]= 0.5*(sedata.s[:,0,0] - sedata.s[:,0,2] + sedata.s[:,2,0] - sedata.s[:,2,2]) #scd11
    scd[:,0,1]= 0.5*(sedata.s[:,0,1] - sedata.s[:,0,3] + sedata.s[:,2,1] - sedata.s[:,2,3]) #scd12
    scd[:,1,0]= 0.5*(sedata.s[:,1,0] - sedata.s[:,1,2] + sedata.s[:,3,0] - sedata.s[:,3,2]) #scd21
    scd[:,1,1]= 0.5*(sedata.s[:,1,1] - sedata.s[:,1,3] + sedata.s[:,3,1] - sedata.s[:,3,3]) #scd22

    return sdd, scc, sdc, scd


# Power Weighting Filter
# A: transmitted peak voltage   
# fb: baudrate
# tr: rise time 7ps
# frequencies defination: freq_start, freq_stop, freq_npoints
def weighting_filter(A, fb, tr, f_start, f_stop, f_npoints):
    freq = np.linspace(f_start, f_stop, f_npoints)
    fr = 0.75*fb
    ft = 0.2635/tr
    # The final weighting fucntion
    w_1 = 1/(1+(freq/ft)**4)
    w_2 = 1/(1+(freq/fr)**8)
    w_3 = (np.sinc(freq/fb))**2 
    w_f = (A**2) * w_1 * w_2 * w_3 
    return w_f, w_1, w_2, w_3, freq


#*********** Main ***********#
# load data using the given SNP file
#sedata = rf.Network('./line1_10in_freq_dep16.s4p')
sedata = rf.Network('./DVN_B1_Multibond_LW_3p6_15in_6.s4p')
#sedata = rf.Network('./DVN_B1_Multibond_LW_3p6_10in_2.s4p')
flag_port1to3_2to4_thru = False # Set the flag to True, if the given s-parameter is 'PORT1 -> PORT3' and 'PORT2 -> PORT4'
# Plot SE insertion loss
sedata.plot_s_db(m=0, n=0)
sedata.plot_s_db(m=1, n=0)
sedata.plot_s_db(m=2, n=0)
sedata.plot_s_db(m=3, n=0)
plt.show()

# modal conversion of the S-parameter
# before the conversion, make sure that the S-parameters's port order
# TX -> RX
# PORT1 -> PORT2
# PORT3 -> PORT4
sdd, scc, sdc, scd = s2smm(sedata, flag_port1to3_2to4_thru)
# Plot differential insertion loss
sdd11 = sdd[:,0,0]
sdd12 = sdd[:,0,1]
sdd21 = sdd[:,1,0]
sdd22 = sdd[:,1,1]
f_start = sedata.frequency.start
f_stop = sedata.frequency.stop
f_npoints = sedata.frequency.npoints
freq = np.linspace(f_start, f_stop, f_npoints)
plt.plot(freq/1e9, 20*np.log10(sdd11), label="SDD11", color='green')
plt.plot(freq/1e9, 20*np.log10(sdd12), label="SDD12", color='red')
plt.plot(freq/1e9, 20*np.log10(sdd21), label="SDD21", color='orange', linestyle='dashed')
plt.plot(freq/1e9, 20*np.log10(sdd22), label="SDD22", color='blue', linestyle='dashed')
plt.title('SDD')
plt.xlabel("Frequency [GHz]")
plt.ylabel("[dB]")
plt.grid(True)
plt.legend()
plt.show()


# Power Weighting Filter
A = 1 # transmitted peak voltage
fb = 32*1e9 #baudrate
#fb = 212.5*1e9 #baudrate
tr = 7*1e-12 #rise time 7ps
#tr = 3*1e-12 #rise time 3ps
w_f, w_1, w_2, w_3, freq = weighting_filter(A, fb, tr, f_start, f_stop, f_npoints)
# plot the weighting filter
plt.plot(freq/1e9, w_1, label="eq1", color='red')
plt.plot(freq/1e9, w_2, label="eq2", color='green')
plt.plot(freq/1e9, w_3, label="eq3", color='blue')
plt.plot(freq/1e9, w_f, label="weighting fucntion", color='orange')
plt.xlabel("Frequency [GHz]")
plt.ylabel("Magnitude")
plt.title("Power Weighting Filter/Components")
plt.grid(True)
plt.legend()
plt.show()

# Retrun loss after the weighting filter
rl_wf = (w_f * (sdd11**2))**0.5
rl_wf_integ = sum(rl_wf) * (1/len(rl_wf))
rl_wf_integ_db = 20*np.log10(np.abs(rl_wf_integ))
print('integrated return loss  = ' + str(round(rl_wf_integ_db,2)) + "dB")
plt.plot(freq/1e9, 20*np.log10(sdd11), label="SDD11", color='green')
plt.plot(freq/1e9, 20*np.log10(rl_wf), label="SDD11_after_WF", color='red', linestyle='dashed')
plt.title('SDD')
plt.xlabel("Frequency [GHz]")
plt.ylabel("[dB]")
plt.grid(True)
plt.legend()
plt.show()
