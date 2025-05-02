import matplotlib.pyplot as plt
import pandas as pd

fig, ax = plt.subplots()
plt.axhline(y=0, color='lightgray', linestyle='-', linewidth=1)

# Color codes for KL bands
# low KL closer to red, high KL closer to violet
color_dict = {1: '#d66808', 2: '#bb4855', 3: '#9e25a9', 4: '#8000ff'}

cols = ['Wavelength (micron)', 'F (flam)', 'Error (flam)']

kl = 50
for day in range(1, 5):
    # read data for current day at given KL
    cur_df = pd.read_csv(f'date{day}/date{day}_spec_KL{kl}.csv', usecols=cols)
    x = cur_df[cols[0]].to_list()
    y = cur_df[cols[1]].to_list()
    err = cur_df[cols[2]].to_list()

    # plot with error
    plt.errorbar(x, y, yerr=err, fmt=' ', ecolor=color_dict[day], elinewidth=0.8)
    plt.plot(x, y, color=color_dict[day])

plt.legend(['_xaxis', '2017-08-31', '2017-09-07', '2018-02-07', '2019-01-13'])

plt.title("All 4 dates at KL 50")
plt.xlabel('$\lambda$ $(\mu m)$')
plt.ylabel('$F_{\lambda}$ $(erg/s/cm^2/\AA)$')
plt.savefig('all_dates_KL50.png')

