from plotting_peak import *

data = read_csv('Data\\15-11-24\dma ab 50um 0min 365nm.Sample.Raw.csv')

plt.plot(data['nm'], data['A'])

plt.show()