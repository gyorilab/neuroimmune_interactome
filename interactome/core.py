import csv
import json
import tqdm
import pandas
import pickle
import obonet
import pystow
from collections import defaultdict
from indra.util import write_unicode_csv
from indra.tools import assemble_corpus as ac
from indra.ontology.bio import bio_ontology
from indra.databases import uniprot_client, hgnc_client
from indra.sources import omnipath
from indra.sources.omnipath.api import _get_interactions
from indra.sources.indra_db_rest import get_curations
from indra.databases import uniprot_client

base_path = pystow.module('neuroimmune')

curated_receptors = {'PLXNA1', 'ITGB2', 'RYR2'}
curated_ligands = {'CRLF1', 'CD274'}
curated_ion_channels = {'P2RX7', 'SCN5A', 'SCN8A', 'TRPC3', 'TRPV6', 'CFTR'}
curated_enzymes = set()


def get_enzymes(cached=True):
    enzymes_path = base_path.join('gene_lists', name='enzymes.csv')
    conflict_curations = get_conflict_curations()
    if not cached or not enzymes_path.exists():
        import pyobo
        ontology = pyobo.get_ontology("eccode")
        obo_path = base_path.join('resources', name='eccode.obo')
        ontology.write_obo(obo_path)
        obo = obonet.read_obo(obo_path)
        up_nodes = {node[8:] for node in obo.nodes if node.startswith('uniprot')}
        human_ups = {u for u in up_nodes if uniprot_client.is_human(u)}
        enzymes = {uniprot_client.get_gene_name(u) for u in human_ups}
        enzymes = {g for g in enzymes if not hgnc_client.is_kinase(g)}
        enzymes = {g for g in enzymes if not hgnc_client.is_phosphatase(g)}
        enzymes = (enzymes | conflict_curations['enzyme']) \
            - curated_receptors - curated_ligands - curated_ion_channels
        with open(base_path.join('gene_lists', name='enzymes.csv'), 'w') as fh:
            fh.write('\n'.join(sorted(enzymes)))
    else:
        with open(enzymes_path) as fh:
            enzymes = fh.read().splitlines()
    return enzymes


def get_omnipath_interactions(cached=True):
    interactions_path = base_path.join('intermediate',
                                       name='omnipath_interactions.json')
    if not cached or not interactions_path.exists():
        ligrec_sources = ['CellPhoneDB', 'Guide2Pharma', 'HPMR', 'ICELLNET',
                          'Kirouac2010', 'CellTalkDB', 'CellChatDB',
                          'connectomeDB2020',
                          'Ramilowski2015', 'talklr']
        interactions = _get_interactions()
        interactions = [i for i in interactions
                        if set(ligrec_sources) & set(i['sources'])]
        interactions = [i for i in interactions
                        if i['curation_effort'] > 0]
        with open(interactions_path, 'w') as fh:
            json.dump(interactions, fh, indent=1)
    else:
        with open(interactions_path, 'r') as fh:
            interactions = json.load(fh)
    return interactions


def get_ligands_receptors(cached=True):
    ligands_path = base_path.join('gene_lists', name='ligands.csv')
    receptors_path = base_path.join('gene_lists', name='receptors.csv')
    cpdb_generated_proteins = base_path.join('resources',
                                             name='protein_generated.csv')
    conflict_curations = get_conflict_curations()
    if not cached or not ligands_path.exists() or not receptors_path.exists():
        interactions = get_omnipath_interactions()
        enzymes = get_enzymes()

        cpdb_df = pandas.read_csv(cpdb_generated_proteins)
        cpdb_receptors_up = cpdb_df[cpdb_df.receptor].uniprot
        cpdb_receptors_hgnc = [uniprot_client.get_gene_name(up)
                               for up in cpdb_receptors_up]
        cpdb_receptors_hgnc = {g for g in cpdb_receptors_hgnc
                               if hgnc_client.get_hgnc_id(g)}
        ligands = ({uniprot_client.get_gene_name(i['source'])
                   for i in interactions if 'COMPLEX' not in i['source']
                   and 'COMPLEX' not in i['target']
                   and i['consensus_direction'] == 1
                   and uniprot_client.is_human(i['source'])} \
            | conflict_curations['ligand']) \
            - curated_receptors - set(enzymes) - curated_ion_channels \
            - curated_ligands - conflict_curations['receptor']
        receptors = (cpdb_receptors_hgnc | \
                  {uniprot_client.get_gene_name(i['target'])
                   for i in interactions if 'COMPLEX' not in i['source']
                   and 'COMPLEX' not in i['target']
                   and i['consensus_direction'] == 1
                   and uniprot_client.is_human(i['target'])} \
            | conflict_curations['receptor']) \
            - curated_ligands - set(enzymes) - curated_ion_channels \
            - conflict_curations['ligand'] - conflict_curations['ion channel'] \
            - conflict_curations['enzyme']

        assert not ligands & receptors

        with open(ligands_path, 'w') as fh:
            fh.write('\n'.join(sorted(ligands)))
        with open(receptors_path, 'w') as fh:
            fh.write('\n'.join(sorted(receptors)))
    else:
        with open(ligands_path) as fh:
            ligands = fh.read().splitlines()
        with open(receptors_path) as fh:
            receptors = fh.read().splitlines()
    return ligands, receptors


def get_ion_channels(cached=True):
    ion_channel_path = base_path.join('gene_lists', name='ion_channels.csv')
    conflict_curations = get_conflict_curations()
    if not cached or not ion_channel_path.exists():
        idg_ion_channels = base_path.join('resources',
                                          name='idg_ion_channels.csv')
        with open(idg_ion_channels) as fh:
            idg_ion_channels = fh.read().splitlines()
        ion_channels = sorted((set(idg_ion_channels) | curated_ion_channels
                               | conflict_curations['ion channel'])
                              - curated_receptors
                              - curated_ligands
                              - curated_enzymes
                              - conflict_curations['receptor'])
        with open(ion_channel_path, 'w') as fh:
            fh.write('\n'.join(ion_channels))
    else:
        with open(ion_channel_path) as fh:
            ion_channels = fh.read().splitlines()
    return ion_channels


def get_ligand_receptor_statements(cached=True):
    ligand_receptor_statements = base_path.join('intermediate',
        name='ligand_receptor_statements.pkl')
    if not cached or not ligand_receptor_statements.exists():
        interactions = get_omnipath_interactions()
        ligands, receptors = get_ligands_receptors()
        op = omnipath.OmniPathProcessor(ligrec_json=interactions)
        op.process_ligrec_interactions()
        statements = [stmt for stmt in op.statements
                      if stmt.agent_list()[0].name in ligands
                      and stmt.agent_list()[1].name in receptors]

        # Enrich Omnipath statements with INDRA evidence
        stmts_by_hash = {stmt.get_hash(): stmt for stmt in statements}
        hashes = list(stmts_by_hash)
        evs_by_hash = download_evidences(hashes)
        for hash, evs in evs_by_hash.items():
            stmts_by_hash[hash].evidence += evs

        dump_stmts_html(statements,
                        base_path.join('intermediate',
                                       name='ligand_receptor_statements.html'))

        with open(ligand_receptor_statements, 'wb') as fh:
            pickle.dump(statements, fh)
    else:
        with open(ligand_receptor_statements, 'rb') as fh:
            statements = pickle.load(fh)
    stmts_to_cpdb(statements, base_path.join('cpdb',
        name='ligand_receptor_interactions.csv'),
        base_path.join('intermediate',
                       name='ligand_receptor_interactions_hgnc.csv'))
    return statements


def get_ligand_ion_channel_statements(cached=True):
    ion_channels = get_ion_channels()
    ligands, _ = get_ligands_receptors()

    indra_sif_path = pystow.join('indra', 'db', name='sif.pkl')
    with open(indra_sif_path, 'rb') as fh:
        indra_df = pickle.load(fh)

    indra_df = indra_df[indra_df['stmt_type'].isin({'Activation', 'Complex'})]
    # Note that we can do this because Complexes show up in both directions
    # in the SIF file so we don't lose any statements
    indra_df = indra_df[indra_df['agA_name'].isin(ligands) &
                        (indra_df['agA_ns'] == 'HGNC') &
                        indra_df['agB_name'].isin(ion_channels) &
                        (indra_df['agB_ns'] == 'HGNC')]

    indra_df = indra_df[indra_df.apply(evidence_filter, axis=1)]

    # Get statements based on hashes from INDRA DB
    hashes = set(indra_df.stmt_hash)
    stmts = list(download_statements(hashes, ev=10000).values())

    # Filter by curation
    curs = get_curations()
    stmts = ac.filter_by_curation(stmts, curs)

    # Filter for directness
    stmts = ac.filter_direct(stmts)
    hashes = {stmt.get_hash() for stmt in stmts}

    indra_df = indra_df[indra_df['stmt_hash'].isin(hashes)]

    indra_df.to_csv(base_path.join('intermediate',
                                   name='ligand_ion_channel_indra_sif.csv'),
                    index=False)

    with open(base_path.join('intermediate',
                             name='ligand_ion_channel_statements.pkl'), 'wb') as fh:
        pickle.dump(stmts, fh)

    dump_stmts_html(stmts,
                    base_path.join('intermediate',
                                   name='ligand_ion_channel_statements.html'))

    rows = [('id_cp_interaction', 'partner_a', 'partner_b', 'source')]
    hgnc_rows = [('id_cp_interaction', 'partner_a', 'partner_b', 'source')]
    interactions = sorted(set(zip(indra_df['agA_name'], indra_df['agB_name'])))
    for idx, (ligand, ion_channel) in enumerate(interactions):
        ion_channel_hgnc = hgnc_client.get_hgnc_id(ion_channel)
        ion_channel_up = hgnc_client.get_uniprot_id(ion_channel_hgnc)
        ligand_hgnc = hgnc_client.get_hgnc_id(ligand)
        ligand_up = hgnc_client.get_uniprot_id(ligand_hgnc)
        if not ligand_up or not ion_channel_up:
            continue
        rows.append(('INDRA-%s' % idx, ligand_up, ion_channel_up, 'INDRA'))
        hgnc_rows.append(('INDRA-%s' % idx, ligand, ion_channel, 'INDRA'))
    write_unicode_csv(
        base_path.join('cpdb', name='ligand_ion_channel_interactions.csv'),
        rows)
    write_unicode_csv(
        base_path.join('intermediate', name='ligand_ion_channel_interactions_hgnc.csv'),
        hgnc_rows)
    return stmts


def get_conflict_curations():
    fname = base_path.join('resources',
                           name='Annotation of conflicted genes_032222.xlsx')
    df = pandas.read_excel(fname, engine='openpyxl', header=None)
    annotations = defaultdict(set)
    for _, row in df.iterrows():
        if not pandas.isna(row[1]) and not pandas.isna(row[2]):
            entity_type = row[2].strip().lower()
            gene_name = row[1].strip()
            annotations[entity_type].add(gene_name)
    return dict(annotations)


def evidence_filter(row):
    readers = {'medscan', 'eidos', 'reach',
               'rlimsp', 'trips', 'sparser', 'isi'}
    if set(row['source_counts']) < readers and row['evidence_count'] == 1:
        return False
    elif set(row['source_counts']) == {'sparser'}:
        return False
    return True


def get_enzyme_product_statements(cached=True):
    pc_sif_path = base_path.join('resources',
                                 name='PathwayCommons12.Detailed.hgnc.sif.gz')
    enzymes = get_enzymes()
    df = pandas.read_csv(pc_sif_path, sep='\t', header=None)
    df = df[df[1] == 'controls-production-of']
    df = df[df[0].isin(set(enzymes))]

    enzymes_by_product = defaultdict(set)
    for _, row in df.iterrows():
        enzyme, product = row[0], row[2]
        chebi_name = bio_ontology.get_name('CHEBI', product)
        enzymes_by_product[chebi_name].add(enzyme)

    products = set(enzymes_by_product) - {'ATP', 'ADP'}

    indra_sif_path = pystow.join('indra', 'db', name='sif.pkl')
    with open(indra_sif_path, 'rb') as fh:
        indra_df = pickle.load(fh)

    _, receptors = get_ligands_receptors()
    ion_channels = get_ion_channels()
    targets = receptors + ion_channels

    indra_df = indra_df[indra_df['stmt_type'].isin({'Activation', 'Complex'})]
    # Note that we can do this because Complexes show up in both directions
    # in the SIF file so we don't lose any statements
    indra_df = indra_df[indra_df['agA_name'].isin(products) &
                        (indra_df['agA_ns'] == 'CHEBI') &
                        indra_df['agB_name'].isin(targets) &
                        (indra_df['agB_ns'] == 'HGNC')]

    indra_df = indra_df[indra_df.apply(evidence_filter, axis=1)]

    # Get statements based on hashes from INDRA DB
    hashes = set(indra_df.stmt_hash)
    stmts = list(download_statements(hashes, ev=10000).values())

    # Filter by curation
    curs = get_curations()
    stmts = ac.filter_by_curation(stmts, curs)

    # Filter for directness
    stmts = ac.filter_direct(stmts)
    hashes = {stmt.get_hash() for stmt in stmts}

    indra_df = indra_df[indra_df['stmt_hash'].isin(hashes)]

    indra_df.to_csv(base_path.join('intermediate',
                                   name='enzyme_product_indra_sif.csv'),
                    index=False)

    with open(base_path.join('intermediate',
                             name='enzyme_product_statements.pkl'), 'wb') as fh:
        pickle.dump(stmts, fh)

    dump_stmts_html(stmts,
        base_path.join('intermediate', name='enzyme_product_statements.html'))

    rows = [('id_cp_interaction', 'partner_a', 'partner_b', 'source')]
    interactions = sorted(set(zip(indra_df['agA_name'], indra_df['agB_name'])))
    idx = 0
    hgnc_rows = [('id_cp_interaction', 'partner_a', 'partner_b', 'source')]
    for product, target in interactions:
        target_hgnc = hgnc_client.get_hgnc_id(target)
        target_up = hgnc_client.get_uniprot_id(target_hgnc)
        enzymes = enzymes_by_product[product]
        for enz in enzymes:
            enz_hgnc = hgnc_client.get_hgnc_id(enz)
            enz_up = hgnc_client.get_uniprot_id(enz_hgnc)
            if not enz_up or not target_up:
                continue
            rows.append(('INDRA-%s' % idx, enz_up, target_up, 'INDRA'))
            hgnc_rows.append(('INDRA-%s' % idx, enz, target, 'INDRA'))
            idx += 1
    write_unicode_csv(base_path.join('cpdb',
                                     name='enzyme_product_interactions.csv'),
                      rows)
    write_unicode_csv(base_path.join('intermediate',
                                     name='enzyme_product_interactions_hgnc.csv'),
                      hgnc_rows)
    return stmts


def stmts_to_cpdb(stmts, cpdb_fname, hgnc_fname):
    rows = [('id_cp_interaction', 'partner_a', 'partner_b', 'source')]
    hgnc_rows = [('id_cp_interaction', 'partner_a', 'partner_b', 'source')]
    added = set()
    for idx, stmt in enumerate(stmts):
        agents = stmt.agent_list()
        up_a = agents[0].db_refs.get('UP')
        up_b = agents[1].db_refs.get('UP')
        if not up_a or not up_b:
            continue
        if (up_a, up_b) in added:
            continue
        rows.append(('INDRA-%s' % idx, up_a, up_b, 'INDRA'))
        hgnc_rows.append(('INDRA-%s' % idx, agents[0].name,
                          agents[1].name, 'INDRA'))
        added.add((up_a, up_b))
    write_unicode_csv(cpdb_fname, rows)
    write_unicode_csv(hgnc_fname, hgnc_rows)


def download_statements(hashes, ev=100):
    """Download the INDRA Statements corresponding to a set of hashes.
    """
    from indra.sources import indra_db_rest
    from indra.util import batch_iter
    stmts_by_hash = {}
    for group in tqdm.tqdm(batch_iter(hashes, 1000),
                           total=int(len(hashes) / 1000)):
        idbp = indra_db_rest.get_statements_by_hash(list(group), ev_limit=ev)
        for stmt in idbp.statements:
            stmts_by_hash[stmt.get_hash()] = stmt
    return stmts_by_hash


def download_evidences(hashes):
    from indra_cogex.client.queries import get_evidences_for_stmt_hashes
    evs_by_hash = get_evidences_for_stmt_hashes(hashes)
    return evs_by_hash


def dump_stmts_html(stmts, fname):
    from indra.assemblers.html import HtmlAssembler
    ha = HtmlAssembler(stmts, db_rest_url='https://db.indra.bio')
    ha.make_model(no_redundancy=True, grouping_level='statement')
    ha.save_model(fname)


def dump_ligand_receptor_omnipath_statements_csv():
    csv_path = base_path.join(name='ligand_receptor_omnipath_statements.csv')
    statements = get_ligand_receptor_statements()
    rows = [('A', 'B', 'Interaction', 'CellPhoneDB', 'Ramilowski2015',
             'sources', 'PMIDs')]
    for stmt in statements:
        pmids = sorted({ev.pmid for ev in stmt.evidence if ev.pmid})
        sources = sorted(stmt.evidence[0].annotations['sources'])
        cpdb = '1' if 'CellPhoneDB' in sources else '0'
        ram = '1' if 'Ramilowski2015' in sources else '0'
        rows.append((stmt.agent_list()[0].name, stmt.agent_list()[1].name,
                     stmt.__class__.__name__, cpdb, ram, ','.join(sources),
                     ','.join(pmids)))
    rows = sorted(rows, key=lambda x: (x[0], x[1], x[5]))

    with open(csv_path, 'w') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerows(rows)


def merge_interactomes():
    interactomes = ['enzyme_product', 'ligand_ion_channel', 'ligand_receptor']
    rows = [['id_cp_interactions', 'partner_a', 'partner_b', 'source']]
    rows_hgnc = [['id_cp_interactions', 'partner_a', 'partner_b', 'source']]
    added = set()
    idx = 0
    for interactome in interactomes:
        this_rows = [['id_cp_interactions', 'partner_a', 'partner_b', 'source']]
        this_rows_hgnc = [['id_cp_interactions', 'partner_a', 'partner_b', 'source']]
        this_idx = 0
        df = pandas.read_csv(base_path.join('cpdb',
                             name='%s_interactions.csv' % interactome))
        df_hgnc = pandas.read_csv(base_path.join('cpdb',
                                  name='%s_interactions_hgnc.csv' % interactome))
        for row, row_hgnc in zip(df.itertuples(), df_hgnc.itertuples()):
            if (row.partner_a, row.partner_b) not in added:
                rows.append(['INDRA-%s' % idx, row.partner_a, row.partner_b,
                             'INDRA'])
                this_rows.append(['INDRA-%s' % this_idx, row.partner_a, row.partner_b,
                                  'INDRA'])
                rows_hgnc.append(['INDRA-%s' % idx, row_hgnc.partner_a,
                                  row_hgnc.partner_b, 'INDRA'])
                this_rows_hgnc.append(['INDRA-%s' % this_idx, row_hgnc.partner_a,
                                      row_hgnc.partner_b, 'INDRA'])
                added.add((row.partner_a, row.partner_b))
                idx += 1
                this_idx += 1
        write_unicode_csv(base_path.join('cpdb', name='%s_unique_interactions.csv' % interactome), this_rows)
        write_unicode_csv(base_path.join('cpdb', name='%s_unique_interactions_hgnc.csv' % interactome),
                          this_rows_hgnc)

    write_unicode_csv(base_path.join('cpdb', name='all_interactions.csv'), rows)
    write_unicode_csv(base_path.join('cpdb', name='all_interactions_hgnc.csv'),
                      rows_hgnc)
