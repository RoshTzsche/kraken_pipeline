"""Run corrected figure analyses, retaining tables and execution status for the paper."""
import argparse
import importlib.util
from pathlib import Path
import pandas as pd
from analysis_support import read_frame, write_run_record
from lefse_support import run_official


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).parent / name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', required=True, help='Full taxonomy matrix, not classification totals')
    parser.add_argument('--metadata', required=True)
    parser.add_argument('--output', default='results/manuscript_review')
    parser.add_argument('--rank', default='genus', help='Primary diversity and ordination rank')
    parser.add_argument('--abundance-rank', default='phylum')
    parser.add_argument('--category', default='Time')
    parser.add_argument('--organism', help='Figure label; defaults to the database name in the input filename')
    parser.add_argument('--threshold', type=float, default=0.01)
    parser.add_argument('--format', choices=['pdf','png','tiff'], default='pdf')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--permutations', type=int, default=999)
    parser.add_argument('--run-lefse', action='store_true', help='Execute installed official LEfSe; otherwise prepare its inputs')
    args = parser.parse_args()
    filename = Path(args.data).stem
    args.organism = args.organism or (filename.removeprefix('Taxonomy_').removesuffix('_Cumulative_Reads')
                                      if filename.startswith('Taxonomy_') else 'Microbiome')
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    frame = read_frame(args.data)
    if 'Rank' not in frame:
        parser.error('The input needs Rank, Name/Scientific Name and sample count columns')
    statuses = []
    def execute(name, function):
        try:
            function()
            statuses.append({'Analysis':name,'Status':'completed','Detail':''})
        except Exception as error:
            statuses.append({'Analysis':name,'Status':'failed','Detail':str(error)})
            print(f'{name}: {error}')
        pd.DataFrame(statuses).to_csv(output/'workflow_status.csv',index=False)
    rare = load_script('07_rarefaction_curve.py')
    bars = load_script('04_generate_Barplots.py')
    alpha = load_script('06_generate_Violin_ANOVA.py')
    beta = load_script('05_generate_PCoA_PieChart.py')
    execute('rarefaction', lambda: rare.generate_rarefaction_curves(args.data,args.rank,None,50,10,
        args.organism,str(output/'01_rarefaction'),args.format,False,args.seed))
    for rank, name in [('domain','02_domains'),(args.abundance_rank,'03_abundance')]:
        if not frame.Rank.str.casefold().eq(rank.casefold()).any():
            statuses.append({'Analysis':name,'Status':'rank_not_available',
                'Detail':f'{rank}: rebuild the matrix from reports using 03_generate_table.py if available'})
            continue
        execute(name, lambda rank=rank,name=name: bars.generate_grouped_microbiome_plots(
            args.data,args.metadata,args.category,rank_level=rank,threshold=args.threshold,
            organism_name=args.organism,output_base=str(output/name),fmt=args.format))
    execute('alpha', lambda: alpha.generate_alpha_diversity_plots(args.data,args.metadata,args.category,
        'SampleID',args.rank,args.organism,str(output/'04_alpha'),args.format,plot_type='boxplot'))
    rank_frame = frame[frame.Rank.str.casefold().eq(args.rank.casefold())].copy()
    execute('beta', lambda: beta.generate_pcoa_plot(rank_frame,args.rank,args.metadata,args.category,'SampleID',
        str(output/'05_pcoa'),args.format,permutations=args.permutations,seed=args.seed,data_path=args.data))
    execute('LEfSe' if args.run_lefse else 'LEfSe_input_preparation',lambda: run_official(rank_frame,
        args.metadata,args.category,'SampleID',str(output/'06_lefse'),args.rank,args.format,
        prepare_only=not args.run_lefse,data_path=args.data))
    pd.DataFrame(statuses).to_csv(output/'workflow_status.csv',index=False)
    write_run_record(output/'workflow','manuscript_figures',[args.data,args.metadata],vars(args))
    if any(row['Status']=='failed' for row in statuses):
        raise SystemExit('Some analyses failed: inspect workflow_status.csv before using figures')
    print(f'Outputs: {output}. LEfSe inputs alone are not an executed LEfSe analysis.')


if __name__=='__main__':
    main()
