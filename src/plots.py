import pandas as pd
import matplotlib.pyplot as plt

methods = ["direct_transfer",
           "meme_on_pssm_searched",
           "meme_on_pssm_searched_seeded",
           "network_transfer",
           "network_transfer_seeded",
           "permuted"]

def read_csv(method):
    return pd.read_csv("../results/%s.csv" % method)

def protein_dist_vs_motif():
    df = read_csv("direct_transfer")
    plt.plot(df['protein_distance'], df['reference_target_dist'], 'o',
                alpha=0.5)
    plt.xlabel("protein distance")
    plt.ylabel("motif-distance")

def method_plot(x, y):
    fig = plt.figure()
    for i in xrange(len(methods)):
        df = read_csv(methods[i])
        ax = fig.add_subplot(2, 3, i+1)
        ax.plot(df[x], df[y], 'o', alpha=0.5)
        ax.set_title(methods[i])
        ax.set_xlabel(x)
        ax.set_ylabel(y)

