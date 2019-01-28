import glob
import numpy as np
import matplotlib as mpl
mpl.use('Qt5Agg')
from matplotlib import pyplot as plt
import pandas
import pandas
from sklearn.manifold import TSNE

colnames = ['theta','intensity']
intensity = []

for file_path in glob.glob('*.csv'):
    with open(file_path,'r') as powder:
        data = pandas.read_csv(powder,names=colnames)
        store = data.theta.tolist()
        iter = len(store)
        intensity.append(store[1:iter])



learn = np.c_[intensity]
print(learn)
learn_embedded = TSNE(n_components=10,perplexity=10).fit_transform(learn)
print(learn_embedded.shape)

plt.scatter(learn_embedded[0],learn_embedded[1])
plt.show()
