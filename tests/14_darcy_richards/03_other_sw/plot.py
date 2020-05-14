import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import colorsys
import yaml

"""
Comparison between Flow123d observe data and S1D code output.
Usage:

  python3 plot.py <out_dir>
  
<out_dir> should contain> flow_observe.yaml, MAIN.IN, W_OBSERV.OUT
"""

# Kelvins law:
# ln(p/p0) = 2 * gamma * M / (rho * r * R * T)
# V_m = M / rho
# M =18
# pressure =  gamma/r ???
#
# 2* gamma * M / r [ N.m/m2 / m * kg/mol ] == (M=18) ???
#
# pressure = -8.314*293.15*998*LN( rh /100) /18
# pressure = -(R = 8.314) *(T=293.15)*(rho=998)*LN( rh /100) /18
#
#pressure = - r_gass_constant * temperature * density * ln( relative_humidity ) / molar_mass

r_gass_constant = 8.314
temperature = 293.15
density = 998
molar_mass = 18
g_const = 9.82

def relative_humidity_element(phead):
    pressure = phead * g_const # 6500 m -> 65000 kPa
    return 100*np.exp( pressure / (r_gass_constant * temperature * density) * molar_mass )


n_obs=8
obs_cols_flow123d=np.array([1,3,4,5,6,7])
#obs_cols_s1d=n_obs-obs_cols_flow123d-1
obs_cols_s1d=obs_cols_flow123d


x = []
output_dir = sys.argv[1]

def s1d_read_obs_positions():
    # read obs positions
    main_file = os.path.join(output_dir, "MAIN.IN")
    with open(main_file, "rt") as f:
        for line in f:
            if "xObs" in line:
                next_line = next(f)
                x_obs_all = list(map(float, next_line.split() ))
                break
    return np.array(x_obs_all)[obs_cols_flow123d[:]]

def s1d_read_obs():       
    x_obs_all = s1d_read_obs_positions()
    
    observ_file = os.path.join(output_dir, "W_OBSERV.OUT")
    with open(observ_file, "rt") as f:    
        raw_table = np.genfromtxt(f, skip_header=1, skip_footer=1)
    times = (raw_table[:,0]).reshape((-1,1)) # keep it matrix
    #sat = raw_table[:,1:n_obs+1]
    
    phead = raw_table[:,11+obs_cols_s1d[:]]
    velocity = -1.0*raw_table[:, 21+obs_cols_s1d[:]]    
    
    return ( np.concatenate( (times, phead), axis=1) , np.concatenate( (times, velocity), axis=1) )   

x_obs= s1d_read_obs_positions()

s1d_data = s1d_read_obs()
phead_s1d = s1d_data[0]
velocity_s1d = s1d_data[1]








# read Flow123d MH results
def read_flow_to_rh(file_name):   
    data_file = os.path.join(output_dir, file_name)
    with open(data_file, "rt") as f:    
        raw_table = np.genfromtxt(f)
    times = (raw_table[:,0]).reshape((-1,1)) # keep it matrix
    #sat = raw_table[:,1:n_obs+1]
    phead = raw_table[:, [2,4,6,8]]
    rh = relative_humidity_element(phead) # drop top and bottom obs.   
    return np.concatenate( (times, rh), axis=1)    

# read Flow123d MH results
def read_flow_yaml_to_rh(file_name):   
    data_file = os.path.join(output_dir, file_name)
    with open(data_file, "rt") as f:
        content=yaml.load(f)
    data = content['data']
    phead=[]
    velocity=[]
    for row in data:        
        time = row['time']
        phead_row = [ time ] + list(np.array(row['pressure_p0'])[obs_cols_flow123d[:]])
        velocity_row = [ time ] + list([ vel[0] for vel in np.array(row['velocity_p0'])[obs_cols_flow123d[:]] ])
        phead.append(phead_row)
        velocity.append(velocity_row)
    phead = np.array(phead)
    #rh = relative_humidity_element(phead[:,1:])
    #rh = np.concatenate( ( phead[:,0:1],  rh ), axis=1)
    velocity = -1.0 * np.array(velocity)[:,1:]
    velocity = np.concatenate( ( phead[:,0:1],  velocity ), axis=1)
    return ( phead, velocity)
        

#rh_flow123d = read_flow_to_rh("result_0.001_corrected_g.txt")
#rh_flow123d = read_flow_to_rh("result_0.001_with_S.txt")
#rh_flow123d=rh_flow123d[:,[0,4,3,2,1]]    

flow_data = read_flow_yaml_to_rh("flow_observe.yaml")

phead_flow123d= flow_data[0]
velocity_flow123d = flow_data[1]

print(x_obs)
print(phead_s1d[1])
print(phead_flow123d[1])


def plot_data(data, title, style, col_hue, point_legend=False):        
    legend_label.append(title)
    legend_handle.append(null_plot)
    for i in range(0,len(x_obs)):
        col_rel = float(i)/len(x_obs)*0.6 + 0.4
        col_rgb=colorsys.hls_to_rgb(col_hue, col_rel, 1.0)        
        ll=plt.plot(data[:,0], data[:,i+1], style)
        plt.setp(ll,linewidth=2, color=col_rgb)
        legend_label.append("")
        legend_handle.append(ll[0])
 
# transpose legend
def transpose_legend(legend, n_cols):
    legend=np.array(legend).reshape(n_cols, -1).transpose()




plt.rcParams.update({'font.size': 12})    
fig = plt.figure(figsize=(15,8))



null_plot=plt.Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
legend_label= [""] + [ str(x_obs[i])+" m" for i in range(0,len(x_obs)) ] 
legend_handle=(len(x_obs)+1)*[null_plot]

plot_data(phead_s1d, "s1d", '--', 0.1)
plot_data(phead_flow123d, "flow", '-',0.1)
#plot_data(relative_humidity_jump, "s1dj", 'm-')
#plot_data(ugn_genuchten, "UGN", '-',0.7)
#plot_data(ugn_brooks, "brooks", 'y-')
#plot_data(ugn_experiment, "exp", '^',1, True)
      
n_col=3
#print np.array(legend_label).reshape(4,-1)
transpose_legend(legend_label, n_col) 
transpose_legend(legend_handle, n_col)

legend = plt.legend(legend_handle, legend_label,  ncol=n_col, handletextpad = -2, frameon=False,
           loc = 'lower right', bbox_to_anchor=(1.0,0.0))

#plt.axis((1000,1500, -200, 10))
plt.axis((00,1500, -2000, 10))
plt.xlabel("time [s]")
plt.ylabel("phead [m]")
#plt.tight_layout(pad=1.08, h_pad=None, w_pad=None, rect=None)
#plt.tight_layout(rec=(0,0.2,1,1))
#plt.tight_layout()
fig.savefig(output_dir + '/compare_phead.svg', bbox_extra_artists=(legend,), bbox_inches='tight')
plt.show()



null_plot=plt.Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
legend_label= [""] + [ str(x_obs[i])+" m" for i in range(0,len(x_obs)) ] 
legend_handle=(len(x_obs)+1)*[null_plot]

plot_data(velocity_flow123d, "flow", '-',0.1)
plot_data(velocity_s1d, "s1d", '--', 0.0)

      
n_col=3
#print np.array(legend_label).reshape(4,-1)
transpose_legend(legend_label, n_col) 
transpose_legend(legend_handle, n_col)

legend = plt.legend(legend_handle, legend_label,  ncol=n_col, handletextpad = -2, frameon=False,
           loc = 'lower right', bbox_to_anchor=(1.0,0.0))

#plt.axis((1000,1500, 0, 1e-6))
plt.axis((0,1500, 0, 1e-4))
plt.xlabel("time [hour]")
plt.ylabel("velocity [m/s]")
#plt.tight_layout(pad=1.08, h_pad=None, w_pad=None, rect=None)
#plt.tight_layout(rec=(0,0.2,1,1))
#plt.tight_layout()
fig.savefig(output_dir + '/compare_velocity.svg', bbox_extra_artists=(legend,), bbox_inches='tight')
plt.show()

