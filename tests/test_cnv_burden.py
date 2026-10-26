import unittest
import numpy as np
from robin.analysis.cnv_burden import analyse_burden as analyse, subtract, packaged_centromeres


class BurdenTests(unittest.TestCase):
    def test_exact_bp_mask_nan_and_chromosome_end(self):
        track = {'chr1': np.array([2., 4., 0., np.nan, 2.])}
        total, _, segments = analyse(track, 100, method='bins', lengths={'chr1': 450}, masks={'chr1': [(125, 175), (150, 190)]})
        self.assertEqual(total['gain_bp'], 35)
        self.assertEqual(total['loss_bp'], 100)
        self.assertEqual(total['neutral_bp'], 150)
        self.assertEqual(total['unknown_bp'], 100)
        self.assertEqual(total['masked_bp'], 65)
        self.assertEqual(total['callable_bp'], 285)
        self.assertAlmostEqual(total['altered_fraction'], 135 / 285)
        self.assertEqual(sum(s['end'] - s['start'] for s in segments), 285)

    def test_difference_neutral_zero_and_negative_loss(self):
        total, _, _ = analyse({'chr1': np.array([0., 1., -1., -2.])}, 10, method='bins', delta=True)
        self.assertEqual((total['gain_bp'], total['loss_bp'], total['neutral_bp']), (10, 20, 10))

    def test_kernel_detects_known_blocks_and_whole_chromosome(self):
        total, _, _ = analyse({'chr1': np.repeat([2., 4., 1.], 100), 'chr2': np.full(100, 4.)}, 10)
        self.assertEqual((total['gain_bp'], total['loss_bp'], total['neutral_bp']), (2000, 1000, 1000))

    def test_genome_totals_use_unbounded_integers(self):
        total, _, _ = analyse({'chr1': np.full(100, 4.), 'chr2': np.full(100, 2.)}, 20_000_000)
        self.assertEqual(total['callable_bp'], 4_000_000_000)
        self.assertIsInstance(total['callable_bp'], int)

    def test_invalid_gap_never_bridged(self):
        total, _, segments = analyse({'chr1': np.array([4., np.nan, 1.])}, 10)
        self.assertEqual(total['unknown_bp'], 10)
        self.assertEqual([(s['start'], s['end']) for s in segments], [(0, 10), (20, 30)])

    def test_no_callable_fails(self):
        with self.assertRaises(ValueError):
            analyse({'chr1': np.array([np.nan, -1.])}, 10)

    def test_robin_default_centromeres(self):
        masks, digest = packaged_centromeres()
        self.assertEqual(len(digest), 64)
        self.assertIn((116796216, 147241828), masks['chr1'])
        track = {'chr1': np.full(200, 2.)}
        totals, _, segments = analyse(track, 1_000_000, method='bins', masks=masks)
        self.assertEqual(totals['masked_bp'], 147241828 - 116796216)
        self.assertTrue(all(s['end'] <= 116796216 or s['start'] >= 147241828 for s in segments))

    def test_overlapping_masks(self):
        self.assertEqual(list(subtract(0, 100, [(10, 30), (20, 40), (50, 200)])), [(0, 10), (40, 50)])


def _sample(tmp_path):
    np.save(tmp_path / 'CNV.npy', {'chr1': np.repeat([2., 4., 1.], 50)})
    np.save(tmp_path / 'CNV_dict.npy', {'bin_width': 1000, 'penalty_value': 10})
    return tmp_path


def test_sources_and_exports(tmp_path):
    import json
    from robin.analysis.cnv_burden import calculate_sample_burden, export_burden
    directory = _sample(tmp_path)
    result = calculate_sample_burden(directory)
    assert result['settings']['source'] == 'sample'
    assert result['summary']['altered_fraction'] == 2/3
    np.save(directory / 'CNV2.npy', {'chr1': np.full(150, 2.)})
    result = calculate_sample_burden(directory)
    assert result['settings']['source'] == 'paired'
    np.save(directory / 'CNV3.npy', {'chr1': np.repeat([0., 2., -1.], 50)})
    result = calculate_sample_burden(directory)
    assert result['settings']['source'] == 'difference'
    assert result['summary']['altered_fraction'] == 2/3
    assert json.loads(export_burden(result))['settings']['centromere_mask_sha256']
    assert b'gain_bp' in export_burden(result, 'summary')
    assert b'median_log2' in export_burden(result, 'segments')


def test_input_replacement_and_deletion_detected(tmp_path):
    from robin.analysis.cnv_burden import input_signature
    directory = _sample(tmp_path)
    before = input_signature(directory)
    (directory / 'CNV.npy').unlink()
    assert before != input_signature(directory)


def test_missing_metadata_and_mismatched_reference_fail(tmp_path):
    import pytest
    from robin.analysis.cnv_burden import calculate_sample_burden
    with pytest.raises(FileNotFoundError):
        calculate_sample_burden(tmp_path)
    directory = _sample(tmp_path)
    np.save(directory / 'CNV2.npy', {'chr1': np.full(149, 2.)})
    with pytest.raises(ValueError, match='must match'):
        calculate_sample_burden(directory, source='paired')


def test_changing_inputs_reject_snapshot(tmp_path, monkeypatch):
    import pytest
    from robin.analysis import cnv_burden as module
    directory = _sample(tmp_path)
    real = module.analyse_burden
    def changing(*args, **kwargs):
        result = real(*args, **kwargs)
        np.save(directory / 'CNV3.npy', {'chr1': np.zeros(150)})
        return result
    monkeypatch.setattr(module, 'analyse_burden', changing)
    with pytest.raises(ValueError, match='changed during calculation'):
        module.calculate_sample_burden(directory)


def test_chart_coordinates_and_fraction(tmp_path):
    from robin.analysis.cnv_burden import calculate_sample_burden
    from robin.gui.components.cnv_burden import burden_chart_options
    result = calculate_sample_burden(_sample(tmp_path))
    profile, bar = burden_chart_options(result)
    assert profile['xAxis']['max'] == .15
    assert profile['series'][0]['connectNulls'] is False
    assert abs(sum(s['data'][0] for s in bar['series']) - 100) < 1e-5
    assert len(profile['series']) == 4


def test_invalid_settings_fail(tmp_path):
    import pytest
    from robin.analysis.cnv_burden import calculate_sample_burden
    directory = _sample(tmp_path)
    for kwargs in ({'gain': -1}, {'loss': 1}, {'penalty': 0}, {'gain': float('nan')}, {'source': 'reference'}):
        with pytest.raises(ValueError):
            calculate_sample_burden(directory, **kwargs)


def test_reference_unmappable_bins_are_unknown(tmp_path):
    from robin.analysis.cnv_burden import calculate_sample_burden
    directory = _sample(tmp_path)
    np.save(directory / 'CNV3.npy', {'chr1': np.zeros(150)})
    reference = np.full(150, 2.)
    reference[10:20] = 0
    np.save(directory / 'CNV2.npy', {'chr1': reference})
    result = calculate_sample_burden(directory)
    assert result['summary']['unknown_bp'] == 10000
    assert result['summary']['callable_bp'] == 140000
    assert result['settings']['reference_mappability_mask'] is True
