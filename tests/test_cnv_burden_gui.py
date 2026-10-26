"""Exercise real NiceGUI controls and the background calculation callback."""
import asyncio
import inspect

import numpy as np
import pytest
from nicegui import Client
from nicegui.page import page

from robin.gui.components.cnv_burden import add_cnv_burden_panel


@pytest.fixture(autouse=True)
def live_timers(monkeypatch):
    import robin.gui.theme as theme
    callbacks = []
    monkeypatch.setattr(theme, "client_timer", lambda interval, callback, **kwargs: callbacks.append(callback))
    return callbacks


def test_panel_calculates_and_invalidates(tmp_path):
    np.save(tmp_path / 'CNV.npy', {'chr1': np.repeat([2., 4., 1.], 50)})
    np.save(tmp_path / 'CNV_dict.npy', {'bin_width': 1000})

    async def exercise():
        with Client(page('/burden-test'), request=None) as client:
            invalidate = add_cnv_burden_panel(tmp_path)
            buttons = [e for e in client.elements.values() if type(e).__name__ == 'Button']
            callback = inspect.getclosurevars(next(iter(buttons[0]._event_listeners.values())).handler).nonlocals['callback']
            assert not buttons[1].enabled
            await callback()
            labels = [e.text for e in client.elements.values() if type(e).__name__ == 'Label']
            assert any('Total altered 66.67%' in text for text in labels)
            assert buttons[1].enabled
            charts = [e for e in client.elements.values() if type(e).__name__ == 'EChart']
            assert all(e.visible for e in charts)
            np.save(tmp_path / 'CNV.npy', {'chr1': np.full(150, 2.)})
            invalidate()
            assert not buttons[1].enabled
            assert not any(e.visible for e in charts)
            await callback()
            labels = [e.text for e in client.elements.values() if type(e).__name__ == 'Label']
            assert any('Total altered 0.00%' in text for text in labels)

    asyncio.run(exercise())


def test_panel_missing_data_reports_error_without_downloads(tmp_path):
    async def exercise():
        with Client(page('/burden-missing'), request=None) as client:
            add_cnv_burden_panel(tmp_path)
            buttons = [e for e in client.elements.values() if type(e).__name__ == 'Button']
            callback = inspect.getclosurevars(next(iter(buttons[0]._event_listeners.values())).handler).nonlocals['callback']
            await callback()
            assert buttons[0].enabled
            assert all(not button.enabled for button in buttons[1:])
            labels = [e.text for e in client.elements.values() if type(e).__name__ == 'Label']
            assert any('Copy Numbr Load unavailable' in text for text in labels)
    asyncio.run(exercise())


def test_cnv_page_contains_burden_panel(tmp_path):
    from types import SimpleNamespace
    from robin.gui.components.cnv import add_cnv_section
    from robin.gui.plotting_preferences import PlottingPreferencesConfig
    with Client(page('/burden-integration'), request=None) as client:
        add_cnv_section(SimpleNamespace(_cnv_state={}, batman_mode=False, plotting_preferences=PlottingPreferencesConfig()), tmp_path)
        labels = [e.text for e in client.elements.values() if type(e).__name__ == 'Label']
        assert 'Copy Numbr Load' in labels
        assert 'Regional CNV events' in labels


def test_live_refresh_updates_overlay_and_metrics(tmp_path, live_timers):
    from nicegui import ui
    np.save(tmp_path / 'CNV.npy', {'chr1': np.repeat([2., 4., 1.], 50)})
    np.save(tmp_path / 'CNV_dict.npy', {'bin_width': 1000})

    async def exercise():
        with Client(page('/load-live'), request=None) as client:
            chart = ui.echart({'series': [{'name': 'CNV', 'type': 'scatter', 'data': [[100, 2]]}]})
            add_cnv_burden_panel(tmp_path, overlay_chart=chart,
                                 plot_context=lambda: ({'chr1': np.zeros(150)}, 1000, 'All'))
            refresh = live_timers[-1]
            await refresh()
            labels = [e.text for e in client.elements.values() if type(e).__name__ == 'Label']
            assert any('Total altered 66.67%' in text for text in labels)
            assert chart.options['series'][0]['name'] == 'CNV'
            areas = [s for s in chart.options['series'] if s['name'] == 'Copy Numbr Load gain'][0]['markArea']['data']
            assert areas[0][0]['xAxis'] == 50000
            assert areas[0][1]['xAxis'] == 100000
            np.save(tmp_path / 'CNV.npy', {'chr1': np.full(150, 2.)})
            await refresh()
            labels = [e.text for e in client.elements.values() if type(e).__name__ == 'Label']
            assert any('Total altered 0.00%' in text for text in labels)
            overlays = [s for s in chart.options['series'] if s['name'].startswith('Copy Numbr Load ')]
            assert len(overlays) == 2
            assert all(not s['markArea']['data'] for s in overlays)
    asyncio.run(exercise())


def test_overlay_offsets_follow_original_plot_and_chromosome_selection():
    from robin.gui.components.cnv_burden import load_overlay_series
    result = {'settings': {'bin_width': 10}, 'segments': [
        {'chrom': 'chr2', 'start': 10, 'end': 30, 'state': 'loss'}]}
    track = {'chr1': np.zeros(5), 'chr2': np.zeros(5), 'chrX': np.zeros(5)}
    series = load_overlay_series(result, track, 10)
    area = series[1]['markArea']['data'][0]
    assert [v['xAxis'] for v in area] == [60, 80]
    area = load_overlay_series(result, track, 10, 'chr2')[1]['markArea']['data'][0]
    assert [v['xAxis'] for v in area] == [10, 30]
    assert not load_overlay_series(result, track, 10, 'chr1')[1]['markArea']['data']
    assert not load_overlay_series(result, track, 20)[1]['markArea']['data']
