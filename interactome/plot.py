import glob
import numpy as np
import pystow
import pandas
import networkx
from random import randint
import matplotlib.pyplot as plt
colors = []

red = '#f2d1da'
green = '#f1fff3'


def get_node_colors(nodes):
    return [green if n in immune_cell_types else red for n in nodes]


def get_node_count(df, node):
    return sum(df['count'][(df['SOURCE'] == node) | (df['TARGET'] == node)])


def normalize_weight(df, w, minw=0.1, maxw=5):
    x = minw + (maxw-minw)*(w - min(df['count'])) / \
        (max(df['count']) - min(df['count']))
    return x


def get_pos(G):
    theta = (np.linspace(0, 1, len(G) + 1)[:-1] * 2 * np.pi).astype(np.float32)
    pos = np.column_stack([np.cos(theta), np.sin(theta)])
    pos = dict(zip(sorted(immune_cell_types) + sorted(neuron_cell_types), pos))
    return pos


base_path = pystow.module('neuroimmune')

if __name__ == '__main__':
    conditions = ['incision', 'UVB', 'zymosan']
    for condition in conditions:
        fname = base_path.join('resources',
                               name='count_network_filtered_%s.txt' % condition)
        df = pandas.read_csv(fname, delimiter='\t')
        immune_cell_types = set(df['SOURCE'])
        neuron_cell_types = set(df['TARGET'])
        G = networkx.from_pandas_edgelist(df,
                                          source='SOURCE',
                                          target='TARGET',
                                          edge_attr='count',
                                          )
        plt.figure(figsize=(11, 11), dpi=300)
        pos = get_pos(G)
        node_labels = {n: n for n in G.nodes()}
        node_labels['Recruited Macs'] = 'Recruited\nMacs'
        p = networkx.draw(G, pos=pos, node_color=get_node_colors(G.nodes()),
                               node_size=[get_node_count(df, node)
                                          for node in G.nodes()],
                               width=[normalize_weight(df, e[2]['count'])
                                      for e in G.edges(data=True)],
                               with_labels=True,
                               edge_color='#4a4d4b',
                               labels=node_labels)
        ax = plt.gca()
        ax.collections[0].set_edgecolor("#bbbbbb")
        plt.subplots_adjust(left=0.225, right=0.9, bottom=0.11, top=0.88)
        plt.savefig(base_path.join('plots', name='network_%s.pdf' % condition))