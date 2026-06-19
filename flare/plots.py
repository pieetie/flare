"""
Save one clean PNG with a strip per metric comparing Flare modes to reference values.
"""
import sys
from .validate import collect

OBS = [('tm', 'Tm  (°C)'), ('self_dimer', 'self dimer  (kcal/mol)'),
       ('cross_dimer', 'cross dimer  (kcal/mol)'), ('hairpin', 'hairpin  (kcal/mol)')]
SERIES = ['flare exact', 'flare literature']
COLORS = {'flare exact': '#2b6cb0', 'flare literature': '#d6453d'}
FOOTER = 'compared against recorded reference values'


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
         + p9.geom_hline(yintercept=0, color='#9aa0a6', size=0.7)
         + p9.geom_jitter(width=0.14, height=0, size=1.7, alpha=0.6,
                          stroke=0, show_legend=True)
         + p9.facet_wrap('observable', scales='free_y', ncol=2)
         + p9.scale_color_manual(values=COLORS, name='')
         + p9.labs(title='flare exact vs flare literature',
                   subtitle='error against reference / one point per oligo / grey line at zero is the reference',
                   caption=FOOTER,
                   x='', y='error  (°C or kcal/mol)')
         + p9.theme_minimal(base_size=12)
         + p9.theme(
             figure_size=(10, 8),
             dpi=150,
             plot_background=p9.element_rect(fill='white', color='white'),
             plot_margin=0.025,
             legend_position='top',
             legend_title=p9.element_blank(),
             legend_box_spacing=0.02,
             legend_key_width=14,
             plot_title=p9.element_text(weight='bold', size=16, ha='center', margin={'b': 4}),
             plot_subtitle=p9.element_text(size=11, color='#5f6368', ha='center', margin={'b': 14}),
             plot_caption=p9.element_text(size=9, color='#80868b', ha='left', margin={'t': 16}),
             strip_text=p9.element_text(weight='bold', size=11, color='#3c4043', margin={'b': 4}),
             axis_title_y=p9.element_text(size=11, color='#5f6368', margin={'r': 10}),
             axis_text_x=p9.element_blank(),
             axis_text_y=p9.element_text(size=9, color='#5f6368'),
             panel_spacing_x=0.10,
             panel_spacing_y=0.09,
             panel_grid_minor=p9.element_blank(),
             panel_grid_major_x=p9.element_blank()))
    g.save(out_png, dpi=150, verbose=False)


def make(path):
    """
    render the error strips straight to the target path
    """
    _chart(_frame(), path)
    return path


def main(argv=None):
    import os
    argv = sys.argv[1:] if argv is None else argv
    path = argv[0] if argv else os.path.join('figures', 'flare_vs_beacon.png')
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    print("wrote", make(path))


if __name__ == "__main__":
    main()
