import numpy as np

f = open('Commands.txt','w')
RodLength = np.linspace(0.9, 1.0, 21)
nv = np.linspace(183, 203, 21);

for i in range(21):
	n = nv[i]
	o = RodLength[i]
	print(o)
	cmdline = 'nohup '+'./simDER' + ' option.txt '  + '--' + ' hangLength ' + '%.4f'% o + ' --' + ' numVertices ' +   '%d' % n +'\n'
	f.write(cmdline)

f.close()
