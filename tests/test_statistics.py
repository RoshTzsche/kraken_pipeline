"""Regression tests for sample inclusion, count denominators and statistical exports."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
sys.path.insert(0, str(SCRIPTS))
from analysis_support import (canonical_sample, count_matrix, metadata_groups,
    compact_letters, display_proportions, rank_comparisons, pcoa_lingoes, dispersion_test)
from lefse_support import prepare_input, parse_results, run_official


def module(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / name)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


rarefaction = module('07_rarefaction_curve.py')
alpha = module('06_generate_Violin_ANOVA.py')
bars = module('04_generate_Barplots.py')
beta = module('05_generate_PCoA_PieChart.py')


class StatisticsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.folder = Path(self.temp.name)
        self.samples = [f's{i}' for i in range(12)]
        rng = np.random.default_rng(123)
        self.frame = pd.DataFrame(rng.integers(1, 200, (5, 12)), columns=self.samples)
        self.frame.insert(0, 'Name', [f'Taxon {i}' for i in range(5)])
        self.frame.insert(0, 'TaxID', range(5))
        self.frame.insert(0, 'Rank', 'Genus')
        self.metadata = self.folder / 'metadata.csv'
        pd.DataFrame({'SampleID': self.samples, 'Time': [0]*4+[7]*4+[30]*4,
                      'Type': ['wild']*4+['lab']*8}).to_csv(self.metadata, index=False)
        self.data = self.folder / 'counts.csv'
        self.frame.to_csv(self.data, index=False)

    def tearDown(self):
        self.temp.cleanup()

    def test_exact_matching_and_confirmed_alias(self):
        self.assertEqual(canonical_sample('NCM2wild_L6_2'), 'ncm2')
        self.assertEqual(canonical_sample('Nwm4_contigs_fasta_proteins'), 'nwm4')
        self.assertNotEqual(canonical_sample('NF10'), canonical_sample('NF1'))

    def test_duplicate_sample_columns_rejected(self):
        with self.assertRaises(ValueError):
            count_matrix(pd.DataFrame({'Rank':['Genus'], 'Ncm2':[1], 'NCM2wild':[2]}))

    def test_missing_counts_are_not_zero(self):
        for value in (np.nan, -1, np.inf):
            with self.assertRaises(ValueError):
                count_matrix(pd.DataFrame({'s1':[value]}))

    def test_unknown_time_and_wild_baseline(self):
        pd.DataFrame({'SampleID':['Ncm2','Nlf4'], 'Time':[np.nan,np.nan],
                      'Type':['wild','lab']}).to_csv(self.metadata,index=False)
        mapping, audit = metadata_groups(['NCM2wild','Nlf4','not_found'],self.metadata)
        self.assertEqual(mapping['NCM2wild'],'0d')
        self.assertIsNone(mapping['Nlf4'])
        self.assertIsNone(mapping['not_found'])
        self.assertEqual(audit.iloc[0].Reason, 'wild_baseline_from_Type')

    def test_conflicting_metadata_rejected(self):
        pd.DataFrame({'SampleID':['s1','S1'], 'Time':[0,7]}).to_csv(self.metadata,index=False)
        with self.assertRaises(ValueError):
            metadata_groups(['s1'],self.metadata)

    def test_other_preserves_denominator_and_missing_samples(self):
        relative = pd.DataFrame({'s1':[.9,.06,.04], 's2':[np.nan]*3}, index=['a','b','c'])
        shown = display_proportions(relative,.1)
        self.assertEqual(shown.loc['a','s1'],90)
        self.assertAlmostEqual(shown.iloc[1,0],10)
        self.assertTrue(shown.s2.isna().all())

    def test_compact_letters_nontransitive_case(self):
        letters=compact_letters(['A','B','C'],[('B','C')])
        self.assertFalse(set(letters['B']) & set(letters['C']))
        self.assertTrue(set(letters['A']) & set(letters['B']))
        self.assertTrue(set(letters['A']) & set(letters['C']))

    def test_rarefaction_seed_and_exact_endpoints(self):
        counts=pd.Series([5,3,2,0])
        one=rarefaction.simulate_rarefaction(counts,[0,1,5,10],10)
        two=rarefaction.simulate_rarefaction(counts,[0,1,5,10],10)
        np.testing.assert_array_equal(one,two)
        np.testing.assert_array_equal(one[0][[0,1,3]],[0,1,3])
        with patch.object(np,'repeat',side_effect=AssertionError('population expansion')):
            rarefaction.simulate_rarefaction(pd.Series([1000000,1]),[100],2)

    def test_rarefaction_expectation(self):
        mean,_=rarefaction.simulate_rarefaction(pd.Series([2,1]),[2],10000)
        self.assertAlmostEqual(mean[0],5/3,delta=.02)

    def test_fractional_rarefaction_counts_rejected(self):
        with self.assertRaises(ValueError):
            rarefaction.simulate_rarefaction(pd.Series([1.5,2]),[1],10)

    def test_alpha_values_and_no_evidence(self):
        result=alpha.compute_alpha_diversity(pd.Series([2,2,0]))
        self.assertAlmostEqual(result['Shannon'],np.log(2))
        self.assertEqual(result['Observed_Richness'],2)
        self.assertTrue(np.isnan(alpha.compute_alpha_diversity(pd.Series([0,0]))['Shannon']))

    def test_no_inference_with_singleton_group(self):
        frame=pd.DataFrame({'Group':['0d','0d','7d'], 'metric':[1,2,3]})
        tests, pairs,_=rank_comparisons(frame,['metric'])
        self.assertEqual(tests.iloc[0].Status,'insufficient_replication')
        self.assertTrue(pairs.empty)

    def test_lingoes_full_coordinates_reconstruct_distances(self):
        d=np.array([[0,1,2,1],[1,0,1,2],[2,1,0,1],[1,2,1,0]],float)
        coords,_,corrected,diagnostic=pcoa_lingoes(d)
        self.assertGreater(diagnostic['lingoes_constant'],0)
        np.testing.assert_allclose(squareform(pdist(coords)),corrected,atol=1e-7)

    def test_permanova_and_permdisp_match_official_implementations(self):
        from skbio import DistanceMatrix
        from skbio.stats.distance import permanova, permdisp
        counts=count_matrix(self.frame).T
        d=squareform(pdist(counts.div(counts.sum(axis=1),axis=0),'braycurtis'))
        groups=['0d']*4+['7d']*4+['30d']*4
        f,r2,p=beta.compute_permanova(d,groups,19,42)
        official=permanova(DistanceMatrix(d),groups,permutations=19,seed=42)
        self.assertEqual(f,official['test statistic'])
        self.assertEqual(p,official['p-value'])
        self.assertTrue(0<=r2<=1)
        dispersion=dispersion_test(d,self.samples,groups,19,42)
        corrected=pcoa_lingoes(d)[2]
        official=permdisp(DistanceMatrix(corrected,ids=self.samples),groups,test='median',permutations=19,seed=42)
        self.assertAlmostEqual(dispersion['Statistic'],official['test statistic'])
        self.assertEqual(dispersion['p_value'],official['p-value'])

    def test_unknown_plot_group_never_enters_inference(self):
        frame=self.frame.copy(); frame['unmatched']=20
        with contextlib.redirect_stdout(io.StringIO()):
            _,tests=beta.generate_pcoa_plot(frame,'genus',str(self.metadata),'Time','SampleID',
                         str(self.folder/'pcoa'),'png',unknown_mode='keep',no_table=True,permutations=19)
        self.assertTrue(tests.N_groups.eq(3).all())
        self.assertTrue(tests.N_samples.eq(12).all())
        run=json.loads((self.folder/'pcoa_run.json').read_text())
        self.assertEqual(run['parameters']['unknown_inference'],'always excluded')

    def test_grouped_bars_average_samples_not_pooled_depths(self):
        pd.DataFrame({'Rank':['Genus','Genus'],'Name':['A','B'],
                      's0':[90,10],'s1':[1,9],'s4':[2,8],'s5':[3,7]}).to_csv(self.data,index=False)
        with contextlib.redirect_stdout(io.StringIO()):
            bars.generate_grouped_microbiome_plots(str(self.data),str(self.metadata),'Time',
                rank_level='genus',threshold=0,output_base=str(self.folder/'bar'),fmt='png',no_table=True)
        shown=pd.read_csv(self.folder/'bar_proportions.csv',index_col=0)
        self.assertAlmostEqual(shown.loc['0d','A'],50)

    def test_lefse_official_input_and_nonhit_parser(self):
        path,audit,labels,groups=prepare_input(self.frame,self.metadata,'Time','SampleID',self.folder/'lefse')
        self.assertEqual(path.read_text().splitlines()[0].split('\t')[1:],groups)
        result=self.folder/'example.res'
        result.write_text('taxid_0\t3.4\t0d\t2.8\t0.003\ntaxid_1\t2.1\t\t\t-\n')
        parsed=parse_results(result,labels)
        self.assertEqual(parsed.iloc[0].LDA_score,2.8)
        self.assertTrue(pd.isna(parsed.iloc[1].LDA_score))

    def test_lefse_missing_dependency_is_explicit(self):
        with patch('lefse_support.shutil.which',return_value=None):
            with self.assertRaisesRegex(RuntimeError,'Official LEfSe'):
                run_official(self.frame,self.metadata,'Time','SampleID',self.folder/'lefse','genus')
        record=json.loads((self.folder/'lefse_run.json').read_text())
        self.assertEqual(record['parameters']['status'],'missing_dependency')

    def test_alpha_and_rarefaction_exports(self):
        with contextlib.redirect_stdout(io.StringIO()):
            alpha.generate_alpha_diversity_plots(str(self.data),str(self.metadata),'Time','SampleID',
                 'genus','Test',str(self.folder/'alpha'),'png',plot_type='boxplot',no_table=True)
            rarefaction.generate_rarefaction_curves(str(self.data),'genus',None,5,3,
                 'Test',str(self.folder/'rarefaction'),'png',True)
        tests=pd.read_csv(self.folder/'alpha_omnibus.csv')
        self.assertEqual(len(tests),4)
        self.assertTrue(tests.p_adjusted.notna().all())
        self.assertTrue((self.folder/'rarefaction_curves.csv').is_file())


if __name__=='__main__':
    unittest.main()
