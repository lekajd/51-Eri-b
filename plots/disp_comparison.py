import matplotlib.pyplot as plt
import pandas as pd

fig, ax = plt.subplots()
plt.axhline(y=0, color='lightgray', linestyle='-', linewidth=1)

# KL values used
kl_values = [1, 5, 10, 15, 20, 25,
             30, 35, 40, 45, 50]

# Color codes for KL bands
# low KL closer to red, high KL closer to violet
kl_colors = {1: 'indianred', 5: 'salmon', 10: 'sandybrown', 15: 'gold',
              20: 'greenyellow', 25: 'mediumspringgreen', 30: 'paleturquoise', 
              35: 'lightskyblue', 40: 'cornflowerblue', 45: 'mediumslateblue', 
              50: 'darkorchid'}

# Color codes for null points
pt_colors = {1:'red', 2:'orange', 3:'gold', 4:'limegreen',
                  5:'dodgerblue', 6:'mediumslateblue', 7:'mediumorchid'}

cols = ['Wavelength (micron)', 'F (flam)']

# graph all null points spectra
for n in range(1, 8):
    for kl in kl_values:
        # read data for that KL band and plot
        cur_df = pd.read_csv(f'data/np{n}/np{n}_spectra_KL{kl}.csv', usecols=cols)
        cur_x = cur_df[cols[0]].tolist()
        cur_y = cur_df[cols[1]].tolist()
        plt.plot(cur_x, cur_y, linewidth=0.5, alpha=0.5, color='gold')

# graph 51 eri b spectra
for kl in kl_values:
    # read data for that KL band and plot
    cur_df = pd.read_csv(f'data/51erib/spec_output_flamspectra_KL{kl}.csv', usecols=cols)
    cur_x = cur_df[cols[0]].tolist()
    cur_y = cur_df[cols[1]].tolist()
    plt.plot(cur_x, cur_y, linewidth=1, color='cornflowerblue')


# axis limits
ax.set_xlim([1.1, 1.8])
ax.set_ylim([-3e-16, 4e-16])

# set up title and labels
plt.title(f"All KL modes for all null points", fontsize=11)
plt.xlabel('$\lambda$ $(\mu m)$')
plt.ylabel('$F_{\lambda}$ $(erg/s/cm^2/\AA)$')
#plt.ylabel('$F_{\lambda}$ $(W/m^2/\mu m)$')
# save and show plot
plt.savefig(f'51eb_vs_npts.png', dpi=300)
plt.show()




