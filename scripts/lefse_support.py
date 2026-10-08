"""Invoke official SegataLab LEfSe; never substitute a custom LDA score."""
from itertools import combinations
from pathlib import Path
import shutil
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from analysis_support import count_matrix, metadata_groups, group_order, write_run_record, read_frame


def prepare_input(frame, metadata, category, sample_id, base, groups_subset=None):
    counts = count_matrix(frame)
    mapping, audit = metadata_groups(counts.columns, metadata, category, sample_id)
    totals = counts.sum(axis=0)
    empty = audit.Sample.isin(totals[totals.eq(0)].index)
    audit.loc[empty, ['Included', 'Reason']] = [False, 'no_taxonomic_detections']
    if groups_subset:
        outside = ~audit.Group.isin(groups_subset) & audit.Included
        audit.loc[outside, ['Included', 'Reason']] = [False, 'outside_requested_comparison']
    samples = audit.loc[audit.Included, 'Sample'].tolist()
    groups = [mapping[s] for s in samples]
    if len(set(groups)) < 2 or min(pd.Series(groups).value_counts()) < 2:
        raise ValueError('LEfSe requires at least two samples in every compared group')
    counts = counts[samples]
    counts = counts.loc[counts.sum(axis=1).gt(0)]
    labels = frame.loc[counts.index, 'Name' if 'Name' in frame else 'Scientific Name'].astype(str)
    features = (['taxid_' + str(x) for x in frame.loc[counts.index, 'TaxID']]
                if 'TaxID' in frame else ['feature_' + str(i) for i in range(len(counts))])
    if len(features) != len(set(features)):
        raise ValueError('Duplicate feature identifiers at the selected rank')
    counts.index = features
    base = Path(base)
    base.parent.mkdir(parents=True, exist_ok=True)
    input_path = Path(str(base) + '_input.tsv')
    with input_path.open('w') as handle:
        handle.write('class\t' + '\t'.join(groups) + '\n')
        counts.to_csv(handle, sep='\t', header=False, index=True)
    pd.DataFrame({'Feature': features, 'Taxon': labels.to_numpy()}).to_csv(str(base) + '_features.csv', index=False)
    counts.to_csv(str(base) + '_counts.csv', index_label='Feature')
    pd.DataFrame({'Sample': samples, 'Group': groups}).to_csv(str(base) + '_groups.csv', index=False)
    return input_path, audit, dict(zip(features, labels)), groups


def parse_results(path, labels):
    rows = []
    for line in Path(path).read_text().splitlines():
        fields = line.split('\t')
        if len(fields) != 5:
            raise ValueError('Unexpected official LEfSe result format: expected five tab-separated fields')
        feature, mean, group, score, p = fields
        rows.append({'Feature': feature, 'Taxon': labels.get(feature, feature),
                     'Log10_Max_Mean': float(mean), 'Enriched_group': group or None,
                     'LDA_score': float(score) if score else np.nan,
                     'KW_p_value': float(p) if p not in ('-', '') else np.nan})
    return pd.DataFrame(rows, columns=['Feature', 'Taxon', 'Log10_Max_Mean',
                                      'Enriched_group', 'LDA_score', 'KW_p_value'])


def run_official(frame, metadata, category, sample_id, base, rank, fmt='pdf',
                 kw_alpha=0.05, wilcox_alpha=0.05, lda_threshold=2.0,
                 top=3, prepare_only=False, groups_subset=None, data_path=None):
    input_path, audit, labels, groups = prepare_input(frame, metadata, category, sample_id, base, groups_subset)
    formatter = shutil.which('lefse_format_input.py')
    runner = shutil.which('lefse_run.py')
    formatted, result_path = str(base) + '.in', str(base) + '.res'
    commands = [[formatter or 'lefse_format_input.py', str(input_path), formatted, '-f', 'r', '-c', '1', '-o', '1000000'],
                [runner or 'lefse_run.py', formatted, result_path, '-a', str(kw_alpha),
                 '-w', str(wilcox_alpha), '-l', str(lda_threshold), '-b', '30', '-f', '0.67', '-y', '1']]
    parameters = {'rank': rank, 'category': category, 'kw_alpha': kw_alpha, 'wilcox_alpha': wilcox_alpha,
                  'lda_threshold': lda_threshold, 'normalisation': 1000000, 'bootstrap_iterations': 30,
                  'bootstrap_fraction': 0.67, 'multiclass_strategy': 'one-against-one',
                  'status': 'prepared_only', 'exploratory': True,
                  'subclasses': 'not specified; no subclass consistency claim',
                  'seed': 'official CLI exposes no seed option; exact stochastic replay not guaranteed'}
    write_run_record(base, 'official_LEfSe', [data_path, metadata], parameters, audit, commands)
    if prepare_only:
        return None
    if not formatter or not runner:
        parameters['status'] = 'missing_dependency'
        write_run_record(base, 'official_LEfSe', [data_path, metadata], parameters, audit, commands)
        raise RuntimeError('Official LEfSe executables not found in PATH. Install SegataLab LEfSe with its R/rpy2 dependencies, or use --prepare-only. Input and provenance were saved.')
    with open(str(base) + '_commands.log', 'w') as log:
        for command in commands:
            proc = subprocess.run(command, text=True, stdout=log, stderr=subprocess.STDOUT, check=False)
            if proc.returncode:
                parameters['status'] = 'failed'
                write_run_record(base, 'official_LEfSe', [data_path, metadata], parameters, audit, commands)
                raise RuntimeError(f'LEfSe exited with {proc.returncode}; inspect {base}_commands.log')
    results = parse_results(result_path, labels)
    results.to_csv(str(base) + '_results.csv', index=False)
    hits = results.dropna(subset=['Enriched_group', 'LDA_score']).sort_values('LDA_score', ascending=False)
    if top > 0:
        hits = hits.groupby('Enriched_group', sort=False).head(top)
    if not hits.empty:
        colors = {g: plt.get_cmap('tab10')(i % 10) for i, g in enumerate(group_order(groups))}
        fig, ax = plt.subplots(figsize=(8, max(3, len(hits) * 0.35)))
        ax.barh(range(len(hits)), hits.LDA_score, color=[colors[g] for g in hits.Enriched_group])
        ax.set_yticks(range(len(hits)), [f'{r.Taxon} ({r.Enriched_group})' for r in hits.itertuples()])
        ax.invert_yaxis()
        ax.set_xlabel('Official LEfSe LDA score')
        ax.set_title(f'Exploratory LEfSe — {rank}')
        fig.tight_layout()
        fig.savefig(str(base) + '.' + ('tif' if fmt == 'tiff' else fmt), dpi=300, bbox_inches='tight')
        plt.close(fig)
    parameters['status'] = 'completed'
    write_run_record(base, 'official_LEfSe', [data_path, metadata], parameters, audit, commands)
    return results


def cli():
    import argparse
    parser = argparse.ArgumentParser(description='Official LEfSe adapter; previous custom scores are superseded.')
    parser.add_argument('-d', '--data', required=True)
    parser.add_argument('-r', '--rank', required=True)
    parser.add_argument('-m', '--metadata', required=True)
    parser.add_argument('-c', '--category', default='Time')
    parser.add_argument('-id', '--sample_id', default='SampleID')
    parser.add_argument('-fmt', '--format', choices=['pdf', 'png', 'tiff'], default='pdf')
    parser.add_argument('-o', '--output')
    parser.add_argument('-org', '--organism', default='Microbiome')
    parser.add_argument('--kw_alpha', type=float, default=0.05)
    parser.add_argument('--wilcox_alpha', type=float, default=0.05)
    parser.add_argument('--lda_threshold', type=float, default=2.0)
    parser.add_argument('--top', type=int, default=3)
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--pairwise', action='store_true')
    parser.add_argument('--cladogram', action='store_true')
    parser.add_argument('--rank_by_lda', action='store_true', help='Compatibility flag: scores always come from official LEfSe')
    parser.add_argument('--label_col', help='Optional taxon label column')
    parser.add_argument('--no_table', action='store_true', help='Evidence CSVs are always retained')
    args = parser.parse_args()
    if args.cladogram:
        parser.error('A rank-only table cannot establish a taxonomy tree. Supply actual lineages to the official cladogram tools.')
    frame = read_frame(args.data)
    frame = frame[frame.Rank.str.casefold().eq(args.rank.casefold())].copy()
    if frame.empty:
        parser.error('Selected rank has no rows')
    if args.label_col:
        if args.label_col not in frame:
            parser.error('Requested label column does not exist')
        frame['Name'] = frame[args.label_col]
    base = Path(args.output or Path(__file__).resolve().parent.parent / 'results' / 'LEfSe' / f'{args.organism}_{args.rank}_LEfSe')
    comparisons = [None]
    if args.pairwise:
        mapping, _ = metadata_groups(count_matrix(frame).columns, args.metadata, args.category, args.sample_id)
        comparisons = list(combinations(group_order(g for g in mapping.values() if g is not None), 2))
    for comparison in comparisons:
        comparison_base = str(base) + ('_' + '_vs_'.join(comparison) if comparison else '')
        run_official(frame, args.metadata, args.category, args.sample_id, comparison_base, args.rank,
                     args.format, args.kw_alpha, args.wilcox_alpha, args.lda_threshold, args.top,
                     args.prepare_only, comparison, args.data)
