"""Live Copy Numbr Load metrics and overlays for the existing CNV plot."""
from __future__ import annotations

import asyncio
import logging
from datetime import datetime
from pathlib import Path

from robin.analysis.cnv_burden import (
    calculate_sample_burden,
    export_burden,
    input_signature,
    packaged_centromeres,
)

COLORS = {"gain": "#c84940", "loss": "#306db1", "neutral": "#aeb9c4"}


def burden_chart_options(result):
    """Return ECharts options; segment display clipping never changes totals."""
    offsets, boundaries, cursor = {}, [], 0
    for row in result["chromosomes"]:
        offsets[row["chrom"]] = cursor
        end = cursor + row["represented_bp"]
        boundaries.append({"xAxis": (cursor + end) / 2e6,
                           "label": {"formatter": row["chrom"][3:], "position": "insideEndTop"}})
        cursor = end
    points = {key: [] for key in COLORS}
    for segment in result["segments"]:
        offset = offsets[segment["chrom"]]
        y = max(-2., min(2., segment["median_log2"]))
        label = f"{segment['chrom']}:{segment['start']:,}–{segment['end']:,}"
        points[segment["state"]].extend([
            {"value": [(offset + segment["start"]) / 1e6, y], "name": label},
            {"value": [(offset + segment["end"]) / 1e6, y], "name": label}, None])
    series = [{"name": key.title(), "type": "line", "data": values,
               "connectNulls": False, "showSymbol": False, "animation": False,
               "lineStyle": {"color": COLORS[key], "width": 2},
               "itemStyle": {"color": COLORS[key]}} for key, values in points.items()]
    masks, _ = packaged_centromeres()
    areas = []
    for row in result["chromosomes"]:
        offset = offsets[row["chrom"]]
        span = row["represented_bp"]
        for start, end in masks.get(row["chrom"], []):
            start, end = max(0, start), min(span, end)
            if start < end:
                areas.append([{"xAxis": (offset+start)/1e6}, {"xAxis": (offset+end)/1e6}])
    series[0]["markArea"] = {"silent": True, "itemStyle": {"color": "rgba(128,128,128,0.16)"}, "data": areas}
    series.append({"name": "Bin signal", "type": "scatter", "symbolSize": 2,
                   "itemStyle": {"color": "#8090a0", "opacity": .35},
                   "data": [[(offsets[chrom]+x)/1e6, y]
                            for chrom, values in result.get("preview", {}).items()
                            for x, y in values], "animation": False, "silent": True})
    settings = result["settings"]
    series[0]["markLine"] = {"silent": True, "symbol": "none", "label": {"show": True},
                             "lineStyle": {"color": "#8895a3", "type": "dashed", "width": .5},
                             "data": boundaries + [{"yAxis": settings["gain_threshold"], "label": {"show": False}},
                                                   {"yAxis": settings["loss_threshold"], "label": {"show": False}}]}
    profile = {"animation": False, "backgroundColor": "transparent",
               "textStyle": {"color": "#8090a0"},
               "legend": {"data": [s["name"] for s in series], "textStyle": {"color": "#8090a0"}},
               "toolbox": {"feature": {"saveAsImage": {"name": "cnv_burden_profile"}, "restore": {}}},
               "tooltip": {"trigger": "item"},
               "grid": {"left": 65, "right": 25, "top": 45, "bottom": 75},
               "xAxis": {"type": "value", "min": 0, "max": cursor/1e6, "name": "Autosomal genomic position (Mb)", "nameLocation": "middle", "nameGap": 32},
               "yAxis": {"type": "value", "min": -2.15, "max": 2.15, "name": "Segment median log2"},
               "dataZoom": [{"type": "inside"}, {"type": "slider", "bottom": 8}], "series": series}
    burden = {"animation": False, "backgroundColor": "transparent", "textStyle": {"color": "#8090a0"},
              "tooltip": {"trigger": "axis", "axisPointer": {"type": "shadow"}},
              "grid": {"left": 25, "right": 25, "top": 15, "bottom": 30},
              "xAxis": {"type": "value", "min": 0, "max": 100, "axisLabel": {"formatter": "{value}%"}},
              "yAxis": {"type": "category", "data": [""], "show": False},
              "series": [{"name": key.title(), "type": "bar", "stack": "burden", "barWidth": 35,
                          "itemStyle": {"color": color},
                          "data": [round(result["summary"][key+"_percent"], 6)]} for key, color in COLORS.items()]}
    return profile, burden


def load_overlay_series(result, cnv_map, bin_width, selected="All"):
    """Gain/loss areas in the existing scatter plot's native bp coordinates."""
    import natsort
    allowed = {f"chr{i}" for i in range(23)} | {"chrX", "chrY"}
    offsets, spans, cursor = {}, {}, 0
    for chrom, values in natsort.natsorted(cnv_map.items()):
        if chrom not in allowed:
            continue
        offsets[chrom] = cursor if selected == "All" else 0
        spans[chrom] = len(values) * bin_width
        cursor += spans[chrom]
    series = []
    for state, color in (("gain", "rgba(200,73,64,0.22)"),
                         ("loss", "rgba(48,109,177,0.22)")):
        areas = []
        if result is not None and result["settings"]["bin_width"] == bin_width:
            for segment in result["segments"]:
                chrom = segment["chrom"]
                if segment["state"] != state or chrom not in offsets:
                    continue
                if selected != "All" and chrom != selected:
                    continue
                start, end = segment["start"], min(segment["end"], spans[chrom])
                if start < end:
                    areas.append([{"name": f"Predicted {state}: {chrom}:{start:,}–{end:,}",
                                   "xAxis": offsets[chrom] + start},
                                  {"xAxis": offsets[chrom] + end}])
        series.append({"name": f"Copy Numbr Load {state}", "type": "scatter", "data": [],
                       "animation": False, "silent": True,
                       "markArea": {"silent": True, "label": {"show": False},
                                    "itemStyle": {"color": color}, "data": areas}})
    return series


def add_cnv_burden_panel(sample_dir, *, overlay_chart=None, plot_context=None):
    """Build live metrics beneath the CNV chart; poll changes off the UI thread."""
    from nicegui import ui
    from robin.classification_config import get_cnv_thresholds
    from robin.gui.theme import client_timer
    from robin.workflow_config import get_cnv_penalty_value

    sample_dir = Path(sample_dir)
    state = {"result": None, "busy": False, "generation": 0, "attempt": None}
    gain_default, loss_default = get_cnv_thresholds("chr1", "XX")
    with ui.column().classes("w-full min-w-0 gap-2"):
        ui.label("Copy Numbr Load").classes("text-lg font-medium")
        ui.label("Gain and loss as a percentage of callable represented autosomal bases. ROBIN centromere/satellite regions are excluded.").classes("classification-insight-meta")
        overlay_enabled = ui.switch("Highlight predicted gains (red) and losses (blue)", value=True)
        with ui.expansion("Load calculation settings").classes("w-full"):
            source = ui.select({"auto": "Auto (prefer CNV3)", "sample": "Sample (CNV)",
                                "difference": "Difference (CNV3)", "paired": "Sample + reference (CNV + CNV2)"}, value="auto", label="Profile").classes("w-64")
            method = ui.select({"kernel": "Changepoint segments", "bins": "Direct bin thresholds"}, value="kernel", label="Method").classes("w-52")
            gain = ui.number("Gain log2 threshold", value=gain_default, step=.05).classes("w-40")
            loss = ui.number("Loss log2 threshold", value=loss_default, step=.05).classes("w-40")
            penalty = ui.number("Changepoint penalty", value=get_cnv_penalty_value(), min=.01).classes("w-40")
        status = ui.label("Waiting for CNV data; load updates automatically.").classes("classification-insight-meta")
        metrics = ui.label("Gain: —   Loss: —   Total altered: —").classes("text-lg font-medium")
        details = ui.label("").classes("classification-insight-meta")
        with ui.expansion("How this is calculated").classes("w-full"):
            ui.label("Linear-kernel changepoints define segments. Their median log2 signal is classified against the gain/loss thresholds. Each unmasked base is counted once. CNV3 is sample minus reference; two copies are added before calculating log2(copy number / 2). Missing or negative copy estimates are excluded. Where CNV2 is available, bins without positive finite reference signal are also excluded. Zero-copy estimates count as potential losses. Thresholds require biological validation; purity and ploidy are not corrected.")
            ui.label("Predicted gain/loss intervals are shaded over the existing copy-number signal and follow its chromosome selection, scale and zoom. Calls use the selected load profile, which may be reference-adjusted. Centromeres and unknown intervals are excluded from load calls. Percentages always describe all callable represented autosomes, even when viewing one chromosome. Without chromosome lengths, terminal bins retain their full width.")
        profile_chart = ui.echart({"series": []}).classes("w-full h-80")
        burden_chart = ui.echart({"series": []}).classes("w-full h-28")
        profile_chart.set_visibility(False)
        burden_chart.set_visibility(False)
        downloads = []

        def apply_overlay():
            if overlay_chart is None or plot_context is None:
                return
            cnv_map, width, selected = plot_context()
            current = overlay_chart.options.get("series", [])
            old = [entry for entry in current
                   if entry.get("name", "").startswith("Copy Numbr Load ")]
            desired = (load_overlay_series(state["result"], cnv_map, width, selected)
                       if overlay_enabled.value and state["result"] is not None else [])
            if old == desired:
                return
            keep = [entry for entry in current
                    if not entry.get("name", "").startswith("Copy Numbr Load ")]
            overlay_chart.options["series"] = keep + desired
            overlay_chart.update()

        def invalidate(message="CNV data or settings changed; updating Copy Numbr Load…"):
            state["generation"] += 1
            state["result"] = None
            state["attempt"] = None
            apply_overlay()
            metrics.set_text("Gain: —   Loss: —   Total altered: —")
            details.set_text("")
            status.set_text(message)
            profile_chart.set_visibility(False)
            burden_chart.set_visibility(False)
            for button in downloads:
                button.disable()

        def check_inputs():
            result = state["result"]
            if result is not None and input_signature(sample_dir) != result["input_signature"]:
                invalidate()
            apply_overlay()

        def download(kind):
            check_inputs()
            result = state["result"]
            if result is not None:
                extension = "json" if kind == "json" else "tsv"
                ui.download.content(export_burden(result, kind), f"copy_numbr_load_{kind}.{extension}")

        async def calculate():
            if state["busy"]:
                return
            invalidate("Updating Copy Numbr Load…")
            state["busy"] = True
            generation = state["generation"]
            state["attempt"] = input_signature(sample_dir)
            calculate_button.disable()
            try:
                result = await asyncio.to_thread(
                    calculate_sample_burden, sample_dir, source=source.value,
                    method=method.value, gain=float(gain.value), loss=float(loss.value),
                    penalty=float(penalty.value))
                if generation != state["generation"] or input_signature(sample_dir) != result["input_signature"]:
                    status.set_text("CNV data changed during calculation; retrying automatically.")
                    return
                options, bar = burden_chart_options(result)
                state["result"] = result
                summary, settings = result["summary"], result["settings"]
                metrics.set_text(f"Gain {summary['gain_percent']:.2f}%   Loss {summary['loss_percent']:.2f}%   Total altered {summary['altered_percent']:.2f}%")
                details.set_text(f"Callable {summary['callable_bp']/1e9:.3f} Gb · Masked {summary['masked_bp']/1e6:.2f} Mb · Unknown {summary['unknown_bp']/1e6:.2f} Mb · {len(result['chromosomes'])}/22 autosomes · {settings['bin_width']:,} bp bins")
                status.set_text(f"Updated {datetime.now().astimezone().strftime('%H:%M:%S %Z')} · {settings['input']} · centromere mask applied.")
                profile_chart.options.clear()
                profile_chart.options.update(options)
                profile_chart.update()
                burden_chart.options.clear()
                burden_chart.options.update(bar)
                burden_chart.update()
                profile_chart.set_visibility(overlay_chart is None)
                apply_overlay()
                burden_chart.set_visibility(True)
                for button in downloads:
                    button.enable()
            except Exception as error:
                logging.exception("Copy Numbr Load calculation failed")
                invalidate(f"Copy Numbr Load unavailable: {error}")
                state["attempt"] = input_signature(sample_dir)
            finally:
                state["busy"] = False
                calculate_button.enable()

        with ui.row().classes("gap-2 flex-wrap"):
            calculate_button = ui.button("Refresh Copy Numbr Load", on_click=calculate)
            for kind, label in (("json", "Download JSON"), ("summary", "Summary TSV"),
                                ("chromosomes", "Chromosomes TSV"), ("segments", "Segments TSV")):
                button = ui.button(label, on_click=lambda kind=kind: download(kind)).props("outline")
                button.disable()
                downloads.append(button)
        async def refresh_live():
            check_inputs()
            if state["busy"]:
                return
            signature = input_signature(sample_dir)
            if state["result"] is None and state["attempt"] != signature:
                # Missing data will be retried when the first file arrives.
                if (sample_dir / "CNV_dict.npy").exists():
                    await calculate()

        overlay_enabled.on_value_change(lambda _: apply_overlay())
        client_timer(2.0, refresh_live)
        for control in (source, method, gain, loss, penalty):
            control.on_value_change(lambda _: invalidate())
    return check_inputs
