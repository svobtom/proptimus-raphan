import json
from glob import glob
from Bio.PDB import PDBParser
import numpy as np
import plotly.graph_objects as go
from plotly.offline import plot
import plotly.io as pio


pdb_parser = PDBParser(QUIET=True)
n_atoms = []
GFNFFca_times = []
PROPTIMUSRAPHANgfnff_times = []


c = 0
for dir in sorted(glob("../calculations/AF*")):
    original_structure = pdb_parser.get_structure("original", glob(f"{dir}/input_PDB/*")[0])
    with open(f'{dir}/comparison.json', 'r', encoding='utf-8') as json_file:
        data = json.load(json_file)
        if data["original / GFN-FFca MAD"] is None: # calculation did not converge
            continue
        GFNFFca_times.append(data["GFN-FFca time"])
        PROPTIMUSRAPHANgfnff_times.append(data["PROPTIMUS RAPHANgfnff time"])
        n_atoms.append(len(list(original_structure.get_residues())))
    c += 1
print(f"{c} files processed.")
print(f"\nAverage speed is {sum(n_atoms)/(sum(PROPTIMUSRAPHANgfnff_times)/3600)} per hour.\n")


fig = go.Figure()
fig.add_trace(go.Scatter(
    x=n_atoms,
    y=GFNFFca_times,
    mode='markers',
    name='GFN-FF<sub>Cα</sub>',
    marker=dict(
        color='red',
        size=10
    )
))

fig.add_trace(go.Scatter(
    x=n_atoms,
    y=PROPTIMUSRAPHANgfnff_times,
    mode='markers',
    name='PROPTIMUS RAPHAN<sub>GFN-FF</sub>',
    marker=dict(
        color='blue',
        size=10
    )
))

x = np.linspace(0, 10000, 100)
y = 1 * x
fig.add_trace(go.Scatter(
    x=x,
    y=y,
    mode='lines',
    name=f"x = y",
    line=dict(
        color='green',
        width=3
    )
))


fig.update_layout(
    width=1200,
    height=1200,
    xaxis=dict(
        title="Number of atoms",
        range=[0, 700],
        title_font=dict(size=35),
        tickfont=dict(size=25),  
        showgrid=True
    ),
    yaxis=dict(
        title="Optimisation time [s]",
        # range=[0, 350000],
        title_font=dict(size=35),
        tickfont=dict(size=25),
        showgrid=True
    ),
    yaxis_type="log",
    legend=dict(
        x=0.01,
        y=0.99,
        traceorder="normal",
        font=dict(size=28),
        bgcolor="rgba(0, 0, 128, 0.1)",
        bordercolor="navy",
        borderwidth=3,
    ),
    hovermode='closest'
)

fig.update_layout(template='simple_white')

plot(fig, filename='plotly_optimization_times.html')
pio.write_image(fig, file='times.pdf')
