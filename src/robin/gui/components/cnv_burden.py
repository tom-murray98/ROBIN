"""On-demand, background CNV burden calculation and interactive GUI plots."""
from __future__ import annotations

import asyncio
import logging
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


def add_cnv_burden_panel(sample_dir):
    """Build panel and return an inexpensive input-invalidation callback."""
    from nicegui import ui
    from robin.classification_config import get_cnv_thresholds
    from robin.workflow_config import get_cnv_penalty_value

    sample_dir = Path(sample_dir)
    state = {"result": None, "busy": False, "generation": 0}
    gain_default, loss_default = get_cnv_thresholds("chr1", "XX")
    with ui.column().classes("w-full min-w-0 gap-2"):
        ui.label("Autosomal CNV burden").classes("text-lg font-medium")
        ui.label("Gain and loss as a percentage of callable represented autosomal bases. ROBIN centromere/satellite regions are excluded.").classes("classification-insight-meta")
        with ui.row().classes("w-full items-center gap-3 flex-wrap"):
            source = ui.select({"auto": "Auto (prefer CNV3)", "sample": "Sample (CNV)",
                                "difference": "Difference (CNV3)", "paired": "Sample + reference (CNV + CNV2)"}, value="auto", label="Profile").classes("w-64")
            method = ui.select({"kernel": "Changepoint segments", "bins": "Direct bin thresholds"}, value="kernel", label="Method").classes("w-52")
            gain = ui.number("Gain log2 threshold", value=gain_default, step=.05).classes("w-40")
            loss = ui.number("Loss log2 threshold", value=loss_default, step=.05).classes("w-40")
            penalty = ui.number("Changepoint penalty", value=get_cnv_penalty_value(), min=.01).classes("w-40")
        status = ui.label("Select Calculate to analyse the current CNV files.").classes("classification-insight-meta")
        metrics = ui.label("Gain: —   Loss: —   Total altered: —").classes("text-lg font-medium")
        details = ui.label("").classes("classification-insight-meta")
        with ui.expansion("How this is calculated").classes("w-full"):
            ui.label("Linear-kernel changepoints define segments. Their median log2 signal is classified against the gain/loss thresholds. Each unmasked base is counted once. CNV3 is sample minus reference; two copies are added before calculating log2(copy number / 2). Missing or negative copy estimates are excluded. Where CNV2 is available, bins without positive finite reference signal are also excluded. Zero-copy estimates count as potential losses. Thresholds require biological validation; purity and ploidy are not corrected.")
            ui.label("The genome plot shows segment medians and sampled bin signal, clipped to ±2 for display; grey shaded intervals are masked. Missing intervals remain gaps. Native-bin resolution is used. Percentages exclude sex chromosomes, and absent chromosomes are not assumed neutral. Without chromosome lengths, terminal bins retain their full width.")
        profile_chart = ui.echart({"series": []}).classes("w-full h-80")
        burden_chart = ui.echart({"series": []}).classes("w-full h-28")
        profile_chart.set_visibility(False)
        burden_chart.set_visibility(False)
        downloads = []

        def invalidate(message="Inputs or settings changed. Calculate again to update burden."):
            state["generation"] += 1
            state["result"] = None
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
                invalidate("CNV files changed. Calculate again to update burden.")

        def download(kind):
            check_inputs()
            result = state["result"]
            if result is not None:
                extension = "json" if kind == "json" else "tsv"
                ui.download.content(export_burden(result, kind), f"cnv_burden_{kind}.{extension}")

        async def calculate():
            if state["busy"]:
                return
            invalidate("Calculating CNV burden…")
            state["busy"] = True
            generation = state["generation"]
            calculate_button.disable()
            try:
                result = await asyncio.to_thread(
                    calculate_sample_burden, sample_dir, source=source.value,
                    method=method.value, gain=float(gain.value), loss=float(loss.value),
                    penalty=float(penalty.value))
                if generation != state["generation"] or input_signature(sample_dir) != result["input_signature"]:
                    status.set_text("Inputs or settings changed during calculation. Calculate again.")
                    return
                options, bar = burden_chart_options(result)
                state["result"] = result
                summary, settings = result["summary"], result["settings"]
                metrics.set_text(f"Gain {summary['gain_percent']:.2f}%   Loss {summary['loss_percent']:.2f}%   Total altered {summary['altered_percent']:.2f}%")
                details.set_text(f"Callable {summary['callable_bp']/1e9:.3f} Gb · Masked {summary['masked_bp']/1e6:.2f} Mb · Unknown {summary['unknown_bp']/1e6:.2f} Mb · {len(result['chromosomes'])}/22 autosomes · {settings['bin_width']:,} bp bins")
                status.set_text(f"Calculated from {settings['input']} ({settings['source']}); centromere mask applied.")
                profile_chart.options.clear()
                profile_chart.options.update(options)
                profile_chart.update()
                burden_chart.options.clear()
                burden_chart.options.update(bar)
                burden_chart.update()
                profile_chart.set_visibility(True)
                burden_chart.set_visibility(True)
                for button in downloads:
                    button.enable()
            except Exception as error:
                logging.exception("CNV burden calculation failed")
                invalidate(f"CNV burden unavailable: {error}")
            finally:
                state["busy"] = False
                calculate_button.enable()

        with ui.row().classes("gap-2 flex-wrap"):
            calculate_button = ui.button("Calculate CNV burden", on_click=calculate)
            for kind, label in (("json", "Download JSON"), ("summary", "Summary TSV"),
                                ("chromosomes", "Chromosomes TSV"), ("segments", "Segments TSV")):
                button = ui.button(label, on_click=lambda kind=kind: download(kind)).props("outline")
                button.disable()
                downloads.append(button)
        for control in (source, method, gain, loss, penalty):
            control.on_value_change(lambda _: invalidate())
    return check_inputs
