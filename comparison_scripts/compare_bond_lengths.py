import math
from glob import glob

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
from numpy.core.defchararray import count
from rdkit import Chem


def statistics(data):
    return f"{round(np.mean(data), 3)}, {round(np.percentile(data, 50), 3)}, {round(np.percentile(data, 75), 3)}, {round(np.percentile(data, 90), 3)}, {round(np.percentile(data, 99), 3)}, {round(np.max(data), 3)}"


def create_histogram(data, title, x_axis_title, y_axis_title, filename):
    fig = go.Figure(data=go.Histogram(x=data,
                                      xbins=dict(start=0,
                                                 size=1)))
    percentil_99 = np.percentile(data, 99)
    fig.add_vline(x=percentil_99,
                  line_width=3,
                  line_dash="dash",
                  line_color="red",
                  annotation_text=f" 99th percentile: {percentil_99:.2f}",
                  annotation_position="top right",
                  annotation_font={"size": 25})
    fig.update_layout(barmode='stack',
                      width=1200,
                      height=700,
                      xaxis_title=x_axis_title,
                      yaxis_title=y_axis_title,
                      title=dict(text=title,
                                 font=dict(size=35),
                                 x=0.5,
                                 xanchor='center'),
                      xaxis = dict(title_font=dict(size=35),
                                   tickfont=dict(size=25),
                                   showgrid=True),
                      yaxis = dict(title_font=dict(size=35),
                                   tickfont=dict(size=25),
                                   showgrid=True),
                      yaxis_type="log",
                      template='simple_white',
                      legend_title_text='Number of<br>side chain<br>hydrogen bonds:',
                      legend=dict(title=dict(font=dict(size=25)),
                                             font=dict(size=25)))
    pio.kaleido.scope.mathjax = None
    pio.write_image(fig, file=filename, engine="kaleido")


def bond_length_differences(pdb_file, ref_pdb_file):
    ref_bond_lengths = {}
    ref_rdkit_structure = Chem.MolFromPDBFile(ref_pdb_file, removeHs=False, sanitize=False)
    ref_conf = ref_rdkit_structure.GetConformer()
    for ref_bond in ref_rdkit_structure.GetBonds():
        i,j = ref_bond.GetBeginAtom().GetIdx(), ref_bond.GetEndAtom().GetIdx()
        if i > j:
            i,j = j,i
        ref_bond_lengths[(i,j)] = math.dist(ref_conf.GetAtomPosition(i), ref_conf.GetAtomPosition(j))

    bond_length_diffs = []
    rdkit_structure = Chem.MolFromPDBFile(pdb_file, removeHs=False, sanitize=False)
    conf = rdkit_structure.GetConformer()
    for bond in rdkit_structure.GetBonds():
        i,j = bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()
        if i > j:
            i,j = j,i
        if (i,j) in ref_bond_lengths.keys():
            bond_length_diffs.append(abs(math.dist(conf.GetAtomPosition(i), conf.GetAtomPosition(j)) - ref_bond_lengths[(i,j)]) * 100) # converted to picometers
    return bond_length_diffs


SETraphan_SETca = []
SEToriginal_SETca = []
SETca_SETraphanca = []
SETraphan_SETraphanca = []
SEToriginal_SETraphanca = []

counter = 0
for calculation_dir in sorted(glob("../calculations/AF-*")):

    try:
        original_structure_file = glob(f"{calculation_dir}/input_PDB/*")[0]
        PROPTIMUS_RAPHANgfnff_structure_file = glob(f"{calculation_dir}/optimised_PDB/*optimised.pdb")[0]
        GFNFFca_structure_file = f"{calculation_dir}/constrained_alpha_carbons_optimisations/original/xtbopt.pdb"
        PROPTIMUS_RAPHANgfnff_GFNFFca_structure_file = f"{calculation_dir}/constrained_alpha_carbons_optimisations/raphan/xtbopt.pdb"

        SETraphan_SETca.extend(bond_length_differences(PROPTIMUS_RAPHANgfnff_structure_file, GFNFFca_structure_file))
        SEToriginal_SETca.extend(bond_length_differences(original_structure_file, GFNFFca_structure_file))
        SETca_SETraphanca.extend(bond_length_differences(GFNFFca_structure_file, PROPTIMUS_RAPHANgfnff_GFNFFca_structure_file))
        SETraphan_SETraphanca.extend(bond_length_differences(PROPTIMUS_RAPHANgfnff_structure_file, PROPTIMUS_RAPHANgfnff_GFNFFca_structure_file))
        SEToriginal_SETraphanca.extend(bond_length_differences(original_structure_file, PROPTIMUS_RAPHANgfnff_GFNFFca_structure_file))

    except OSError:
        continue
    counter += 1


print("Comparison of bond lengths:\n")
print(f"{counter} structures processed")
print(f"{len(SEToriginal_SETraphanca)} bonds processed\n")
print("                                       Mean, 50th percentile, 75th percentile, 90th percentile, 99th percentile, maximum value")
print(f"Statistics(SETraphan,          SETca): {statistics(SETraphan_SETca)}")
print(f"Statistics(SEToriginal,        SETca): {statistics(SEToriginal_SETca)}")
print("--------------------------------------------------------------------------------")
print(f"Statistics(SETraphan,   SETraphan+ca): {statistics(SETraphan_SETraphanca)}")
print(f"Statistics(SETca,       SETraphan+ca): {statistics(SETca_SETraphanca)}")
print(f"Statistics(SEToriginal, SETraphan+ca): {statistics(SEToriginal_SETraphanca)}")


create_histogram(SETraphan_SETca, "Bond length deviations between SET<sub>RAPHAN</sub> and SET<sub>Cα</sub>", 'Deviation of bond lengths [pm]', 'Number of bonds', "SETraphan_SETca_bond_lengths.pdf")
create_histogram(SEToriginal_SETca, "Bond length deviations between SET<sub>ORIG</sub> and SET<sub>Cα</sub>", 'Deviation of bond lengths [pm]', 'Number of bonds', "SEToriginal_SETca_bond_lengths.pdf")
create_histogram(SETca_SETraphanca, "Bond length deviations between SET<sub>Cα</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of bond lengths [pm]', 'Number of bonds', "SETca_SETraphan+ca_bond_lengths.pdf")
create_histogram(SETraphan_SETraphanca, "Bond length deviations between SET<sub>RAPHAN</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of bond lengths [pm]', 'Number of bonds', "SETraphan_SETraphan+ca_bond_lengths.pdf")
create_histogram(SEToriginal_SETraphanca, "Bond length deviations between SET<sub>ORIG</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of bond lengths [pm]', 'Number of bonds', "SEToriginal_SETraphan+ca_bond_lengths.pdf")
