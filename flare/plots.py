"""
Save one clean PNG with a curve per metric comparing Flare modes to reference values.
"""
import os
import sys
import tempfile
from .validate import collect, stats

OBS = [('tm', 'Tm  (C)'), ('self_dimer', 'self dimer  (kcal/mol)'),
       ('cross_dimer', 'cross dimer  (kcal/mol)'), ('hairpin', 'hairpin  (kcal/mol)')]
SERIES = ['flare exact', 'flare literature']
COLORS = {'flare exact': '#2b6cb0', 'flare literature': '#d6453d'}


def _frame():
    """
    One row per oligo and Flare mode with the error versus the reference value.
    """
    import pandas as pd
    lit, exact = collect('literature'), collect('beacon_exact')
    rec = []
    for obs, olabel in OBS:
        for rl, rb in zip(lit[obs], exact[obs]):
            rec.append({'observable': olabel, 'series': 'flare exact', 'error': round(rb[-1] - rb[-2], 3)})
            rec.append({'observable': olabel, 'series': 'flare literature', 'error': round(rl[-1] - rl[-2], 3)})
    df = pd.DataFrame(rec)
    df['observable'] = pd.Categorical(df['observable'], [o[1] for o in OBS])
    df['series'] = pd.Categorical(df['series'], SERIES)
    return df


def _chart(df, out_png):
    """
    Two strips of error per metric, one per Flare mode; zero is the reference value.
    """
    import plotnine as p9
    g = (p9.ggplot(df, p9.aes('series', 'error', color='series'))
         + p9.geom_hline(yintercept=0, color='#888888', size=0.7)
         + p9.geom_jitter(width=0.16, height=0, size=1.5, alpha=0.55)
         + p9.facet_wrap('observable', scales='free_y', ncol=2)
         + p9.scale_color_manual(values=COLORS, name='')
         + p9.labs(title='flare exact vs flare literature',
                   subtitle='error against reference, one point per oligo, the grey line at zero is reference',
                   x='', y='error  (C or kcal/mol)')
         + p9.theme_minimal(base_size=12)
         + p9.theme(figure_size=(9.5, 7.5), legend_position='top',
                    plot_title=p9.element_text(weight='bold', size=14),
                    axis_text_x=p9.element_blank(),
                    panel_spacing=0.04))
    g.save(out_png, dpi=150, verbose=False)


def _add_source_footer(chart_png, path):
    """
    Add a bottom line noting that values are compared to recorded references.
    """
    from PIL import Image, ImageDraw, ImageFont
    from matplotlib import font_manager
    chart = Image.open(chart_png).convert('RGBA')
    band = 40
    out = Image.new('RGBA', (chart.width, chart.height + band), 'white')
    out.paste(chart, (0, 0), chart)
    draw = ImageDraw.Draw(out)
    font = ImageFont.truetype(font_manager.findfont('DejaVu Sans'), 15)
    y = chart.height + band // 2
    draw.text((14, y), "compared against recorded reference values",
              fill='#555555', font=font, anchor='lm')
    out.convert('RGB').save(path)


def make(path):
    """
    render the curves then drop the source footer on the image
    """
    with tempfile.TemporaryDirectory() as tmp:
        chart = os.path.join(tmp, 'chart.png')
        _chart(_frame(), chart)
        _add_source_footer(chart, path)
    return path


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    path = argv[0] if argv else "flare_validation.png"
    print("wrote", make(path))


if __name__ == "__main__":
    main()
