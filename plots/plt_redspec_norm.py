import matplotlib.pyplot as plt
import pandas as pd

fig, ax = plt.subplots()
plt.axhline(y=0, color='lightgray', linestyle='-', linewidth=1)

# Color codes for KL bands
# low KL closer to red, high KL closer to violet
color_dict = {'H': 'silver', 'IFS_YH': 'plum', 'IFS_YJ': 'lightblue'}

published_cols = ['wavelength', 'F_lam', 'uncertainty']

# CHARIS data

charis_cols = ['Wavelength (micron)', 'F (flam)']
charis_df = pd.read_csv(f'charis_data/spec_output_flamspectra_KL50.csv', usecols=charis_cols)
x = charis_df[charis_cols[0]].to_list()

redspec_cols = ['wvs (microns)', 'wvs (angstroms)', 'Flam (ergs/cm^2/s/AA)', 'Error (ergs/cm^2/s/AA)']
redspec_df = pd.read_excel('redspec_converted.xlsx', usecols=redspec_cols)
y = redspec_df[redspec_cols[2]].to_list()
y_converted = [(value * 10) for value in y]
err = redspec_df[redspec_cols[3]].to_list()
err_converted = [(value * 10) for value in err]

# normalization & scaling
norm_factor = y_converted[3]
y_normalized = [(value / norm_factor) for value in y_converted]
err_normalized = [(value / norm_factor) for value in err_converted]

scale_factor = 8.301130009117527e-17
y_scaled = [(value * scale_factor) for value in y_normalized]
err_scaled = [(value * scale_factor) for value in err_normalized]

print(x)
print(y_scaled)
print(err_scaled)

#y = charis_df[charis_cols[1]].to_list()
#y_converted = [(value * 10) for value in y]
plt.errorbar(x, y_scaled, yerr=err_scaled, linewidth=1.5, color='dodgerblue', zorder=10,
             elinewidth=1, ecolor='dodgerblue', capsize=2, capthick=1)


## GPI & SPHERE data

datasets = ['H', 'IFS_YH', 'IFS_YJ']
for cur in datasets:
    # read data for current dataset
    cur_df = pd.read_csv(f'51erib_{cur}.csv', usecols=published_cols)
    x = cur_df[published_cols[0]].to_list()
    y = cur_df[published_cols[1]].to_list()
    err = cur_df[published_cols[2]].to_list()

    # plot with error
    plt.errorbar(x, y, yerr=err, fmt=' ', marker='s', ms=5, color=color_dict[cur],
                 capsize=2, capthick=1, ecolor=color_dict[cur], elinewidth=1)


plt.legend(['_xaxis', 'CHARIS (norm.)', 'GPI H', 'SPHERE YH', 'SPHERE YJ'])

plt.title("CHARIS (norm.) & published data")
plt.xlabel('$\lambda$ $(\mu m)$')
plt.ylabel('$F_{\lambda}$ $(W/m^2/\mu m)$')
plt.savefig('redspec.png')



