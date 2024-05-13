import os
import pickle
import pystow
import pandas
from matplotlib_venn import venn3
import matplotlib.pyplot as plt


path = '/Users/ben/hmsdropbox/singlecell_model_analysis/panacea_indra/' \
       'cellphonedb_processor/cpdb_output/all_interactome_AJ analysis'


models = ['incision', 'uvb', 'zymosan']
# This is a file that contains columns to keep from the interaction files
# and serves as a filter over cell types to consider interactions between
# (to avoid neuron-neuron interactions, etc.) An earlier version of this
# file contains column names in the first row cols to keep.xlsx, the more
# recent cols_to_keep_all.xlsx has column names in the first column but with
# gaps
cols_to_keep_sheet = os.path.join(path, 'cols_to_keep_all.xlsx')

base_path = pystow.module('neuroimmune')

if __name__ == '__main__':
    df = pandas.read_excel(cols_to_keep_sheet, sheet_name='Sheet1', header=None)
    keep_columns = list(df[0][~pandas.isna(df[0])])

    interaction_cols = [c for c in keep_columns if '|' in c]

    pkls = ['ligand_receptor_statements.pkl',
            'ligand_ion_channel_statements.pkl',
            'enzyme_product_statements.pkl']
    stmts = []
    for pkl in pkls:
        with open(base_path.join('intermediate', name=pkl), 'rb') as fh:
            stmts += pickle.load(fh)

    indra_stmts = [s for s in stmts if any(ev.source_api != 'omnipath'
                                           for ev in s.evidence)]
    indra_only_stmts = [s for s in indra_stmts
                        if all(ev.source_api != 'omnipath'
                               for ev in s.evidence)]

    sig_ints = {}
    indra_ints = {}
    indra_only_ints = {}
    for model in models:
        fname = os.path.join(path, model, model,
                             'significant_means_%s.txt' % model)
        df = pandas.read_csv(fname, sep='\t')
        df = df[keep_columns]
        df['sum'] = df[interaction_cols].sum(axis=1)

        df_sig = df[df['sum'] > 0]
        sig_ints[model] = set(df_sig['interacting_pair'])

        indra_ints[model] = [
            s for s in indra_stmts if
            (('_'.join([a.name for a in s.agent_list()]) in sig_ints[model])
              or
             ('_'.join([a.name for a in [s.agent_list()[1], s.agent_list()[0]]]))
             in sig_ints[model])
            ]
        indra_only_ints[model] = [
            s for s in indra_only_stmts if
            (('_'.join([a.name for a in s.agent_list()]) in sig_ints[model])
             or
             ('_'.join([a.name for a in [s.agent_list()[1], s.agent_list()[0]]]))
             in sig_ints[model])
        ]
    plt.figure()
    names = ['Incision', 'UV burn', 'Zymosan']
    colors = ['r', 'b', 'g']
    venn3([sig_ints[model] for model in models],
          set_labels=names,
          set_colors=colors)
    plt.savefig('venn_diagram_v2.pdf')
