from glob import glob

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
from rdkit import Chem
from rdkit.Chem import rdMolTransforms


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


def get_dihedral_angles(pdb_file, ref_pdb_file):
    ref_angles = {}
    ref_rdkit_structure = Chem.MolFromPDBFile(ref_pdb_file, removeHs=False, sanitize=False)
    ref_conf = ref_rdkit_structure.GetConformer()
    for xtb_bond in ref_rdkit_structure.GetBonds():
        i,j = xtb_bond.GetBeginAtom().GetIdx(), xtb_bond.GetEndAtom().GetIdx()
        if i > j:
            i,j = j,i
        i_bonded_atoms = [a.GetIdx() for a in ref_rdkit_structure.GetAtomWithIdx(i).GetNeighbors()]
        j_bonded_atoms = [a.GetIdx() for a in ref_rdkit_structure.GetAtomWithIdx(j).GetNeighbors()]
        if len(i_bonded_atoms) and len(j_bonded_atoms):
            for ia in i_bonded_atoms:
                for ja in j_bonded_atoms:
                    if len(set((ia,i,j,ja))) != 4:
                        continue
                    ref_angles[((ia,i,j,ja))] = rdMolTransforms.GetDihedralDeg(ref_conf, ia,i,j,ja)

    angles_diffs = []
    rdkit_structure = Chem.MolFromPDBFile(pdb_file, removeHs=False, sanitize=False)
    conf = rdkit_structure.GetConformer()
    for xtb_bond in rdkit_structure.GetBonds():
        i,j = xtb_bond.GetBeginAtom().GetIdx(), xtb_bond.GetEndAtom().GetIdx()
        if i > j:
            i,j = j,i
        i_bonded_atoms = [a.GetIdx() for a in rdkit_structure.GetAtomWithIdx(i).GetNeighbors()]
        j_bonded_atoms = [a.GetIdx() for a in rdkit_structure.GetAtomWithIdx(j).GetNeighbors()]
        if len(i_bonded_atoms) and len(j_bonded_atoms):
            for ia in i_bonded_atoms:
                for ja in j_bonded_atoms:
                    if len(set((ia,i,j,ja))) != 4:
                        continue
                    if (ia,i,j,ja) in ref_angles.keys():
                        angles_diffs.append(abs(rdMolTransforms.GetDihedralDeg(conf, ia,i,j,ja) - ref_angles[(ia,i,j,ja)]))
    return angles_diffs


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
    
        SETraphan_SETca.extend(get_dihedral_angles(PROPTIMUS_RAPHANgfnff_structure_file, GFNFFca_structure_file))
        SEToriginal_SETca.extend(get_dihedral_angles(original_structure_file, GFNFFca_structure_file))
        SETca_SETraphanca.extend(get_dihedral_angles(GFNFFca_structure_file, PROPTIMUS_RAPHANgfnff_GFNFFca_structure_file))
        SETraphan_SETraphanca.extend(get_dihedral_angles(PROPTIMUS_RAPHANgfnff_structure_file, PROPTIMUS_RAPHANgfnff_GFNFFca_structure_file))
        SEToriginal_SETraphanca.extend(get_dihedral_angles(original_structure_file, PROPTIMUS_RAPHANgfnff_GFNFFca_structure_file))

    except OSError:
        continue
    counter += 1

SETraphan_SETca_transformed = []
SEToriginal_SETca_transformed = []
SETca_SETraphanca_transformed = []
SETraphan_SETraphanca_transformed = []
SEToriginal_SETraphanca_transformed = []

for value in SETraphan_SETca:
    transformed_value = value
    if value > 180:
        transformed_value = 360 - value
    SETraphan_SETca_transformed.append(transformed_value)

for value in SEToriginal_SETca:
    transformed_value = value
    if value > 180:
        transformed_value = 360 - value
    SEToriginal_SETca_transformed.append(transformed_value)

for value in SETca_SETraphanca:
    transformed_value = value
    if value > 180:
        transformed_value = 360 - value
    SETca_SETraphanca_transformed.append(transformed_value)

for value in SETraphan_SETraphanca:
    transformed_value = value
    if value > 180:
        transformed_value = 360 - value
    SETraphan_SETraphanca_transformed.append(transformed_value)

for value in SEToriginal_SETraphanca:
    transformed_value = value
    if value > 180:
        transformed_value = 360 - value
    SEToriginal_SETraphanca_transformed.append(transformed_value)

print("Comparison of dihedral angles:\n")
print(f"{counter} structures processed")
print(f"{len(SEToriginal_SETraphanca)} dihedral angles processed\n")
print("                                       Mean, 50th percentile, 75th percentile, 90th percentile, 99th percentile, maximum value")
print(f"Statistics(SETraphan,          SETca): {statistics(SETraphan_SETca)}")
print(f"Statistics(SEToriginal,        SETca): {statistics(SEToriginal_SETca)}")
print("--------------------------------------------------------------------------------")
print(f"Statistics(SETraphan,   SETraphan+ca): {statistics(SETraphan_SETraphanca)}")
print(f"Statistics(SETca,       SETraphan+ca): {statistics(SETca_SETraphanca)}")
print(f"Statistics(SEToriginal, SETraphan+ca): {statistics(SEToriginal_SETraphanca)}")

create_histogram(SETraphan_SETca_transformed, "Dihedral angle deviations between SET<sub>RAPHAN</sub> and SET<sub>Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SETraphan_SETca_dihedral_angles.pdf")
create_histogram(SEToriginal_SETca_transformed, "Dihedral angle deviations between SET<sub>ORIG</sub> and SET<sub>Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SEToriginal_SETca_dihedral_angles.pdf")
create_histogram(SETca_SETraphanca_transformed, "Dihedral angle deviations between SET<sub>Cα</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SETca_SETraphan+ca_dihedral_angles.pdf")
create_histogram(SETraphan_SETraphanca_transformed, "Dihedral angle deviations between SET<sub>RAPHAN</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SETraphan_SETraphan+ca_dihedral_angles.pdf")
create_histogram(SEToriginal_SETraphanca_transformed, "Dihedral angle deviations between SET<sub>ORIG</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SEToriginal_SETraphan+ca_dihedral_angles.pdf")

