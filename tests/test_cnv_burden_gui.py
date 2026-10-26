"""Exercise real NiceGUI controls and the background calculation callback."""
import asyncio
import inspect

import numpy as np
from nicegui import Client
from nicegui.page import page

from robin.gui.components.cnv_burden import add_cnv_burden_panel


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
            assert any('CNV burden unavailable' in text for text in labels)
    asyncio.run(exercise())


def test_cnv_page_contains_burden_panel(tmp_path):
    from types import SimpleNamespace
    from robin.gui.components.cnv import add_cnv_section
    from robin.gui.plotting_preferences import PlottingPreferencesConfig
    with Client(page('/burden-integration'), request=None) as client:
        add_cnv_section(SimpleNamespace(_cnv_state={}, batman_mode=False, plotting_preferences=PlottingPreferencesConfig()), tmp_path)
        labels = [e.text for e in client.elements.values() if type(e).__name__ == 'Label']
        assert 'Autosomal CNV burden' in labels
        assert 'Regional CNV events' in labels
