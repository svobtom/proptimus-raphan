import math
from collections import defaultdict
from glob import glob

import biotite.structure as struc
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
from Bio.PDB import PDBParser, Superimposer
from biotite.structure import io as biotite


def statistics(data):
    return f"{round(np.mean(data), 3)}, {round(np.percentile(data, 50), 3)}, {round(np.percentile(data, 75), 3)}, {round(np.percentile(data, 90), 3)}, {round(np.percentile(data, 99), 3)}, {round(np.max(data), 3)}"


def create_histogram_hbonds(structures, diff, diffs, file, title=None):
    num_of_H_bonds_diffs = defaultdict(list)
    for structure in structures:
        for res in structure.get_residues():
            num_of_side_chain_H_bonds = len([atom for atom in res.get_atoms() if atom._h_bond])
            num_of_H_bonds_diffs[num_of_side_chain_H_bonds].extend([getattr(atom, diff) for atom in res.get_atoms()])

    data = []
    for color, (key, value) in zip(["#AB63FA", "#636EFA", "#19D3F3", "#00CC96", "#B6E880", "#FECB52", "#FFA15A", "#EF553B"], sorted(num_of_H_bonds_diffs.items(), reverse=True)):
        data.append(go.Histogram(x=value,
                                 name=key,
                                 marker=dict(color=color),
                                 xbins=dict(start=0,
                                            size=0.1)))
    fig = go.Figure(data=data)
    percentil_99 = np.percentile(diffs, 99)
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
                      xaxis_title='Deviation of atomic positions [Å]',
                      yaxis_title='Number of atoms',
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
    if title:
        fig.update_layout(title=dict(text=title,
                                     font=dict(size=35),
                                     x=0.5,
                                     xanchor='center'))
    else:
        file = file.replace(".pdf", "_without_title.pdf")
    pio.kaleido.scope.mathjax = None
    pio.write_image(fig, file=file, engine="kaleido")

def create_histogram(data, title, x_axis_title, y_axis_title, filename):
    fig = go.Figure(data=go.Histogram(x=data,
                                      xbins=dict(start=0,
                                                 size=0.1)))
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


def atom_coordinates_diferences(prot1, prot2):
    sup.set_atoms([a for a in prot1.get_atoms() if a.name == "CA"], [a for a in prot2.get_atoms() if a.name == "CA"])
    sup.apply(prot2.get_atoms())
    d = []
    for a1, a2 in zip(prot1.get_atoms(), prot2.get_atoms()):
        if a1.full_id != a2.full_id:
            exit("ERROR! Compared atoms are not identical!")
        d.append(math.dist(a1.coord, a2.coord))
    return d


SETraphan_SETca = []
SEToriginal_SETca = []
SETca_SETraphanca = []
SETraphan_SETraphanca = []
SEToriginal_SETraphanca = []

sup = Superimposer()
pdb_parser = PDBParser(QUIET=True)
counter = 0
structures = []

for calculation_dir in sorted(glob("../calculations/AF-*")):
    try:
        original_structure = pdb_parser.get_structure("structure", glob(f"{calculation_dir}/input_PDB/*")[0])
        PROPTIMUS_RAPHANgfnff_structure = pdb_parser.get_structure("structure", glob(f"{calculation_dir}/optimised_PDB/*optimised.pdb")[0])
        GFNFFca_structure = pdb_parser.get_structure("structure", f"{calculation_dir}/constrained_alpha_carbons_optimisations/original/xtbopt.pdb")
        PROPTIMUS_RAPHANgfnff_GFNFFca_structure = pdb_parser.get_structure("structure", f"{calculation_dir}/constrained_alpha_carbons_optimisations/raphan/xtbopt.pdb")

        SETraphan_SETca_diffs = atom_coordinates_diferences(GFNFFca_structure, PROPTIMUS_RAPHANgfnff_structure)
        SEToriginal_SETca_diffs = atom_coordinates_diferences(GFNFFca_structure, original_structure)
        SETca_SETraphanca_diffs = atom_coordinates_diferences(GFNFFca_structure, PROPTIMUS_RAPHANgfnff_GFNFFca_structure)
        SETraphan_SETraphanca_diffs = atom_coordinates_diferences(PROPTIMUS_RAPHANgfnff_structure, PROPTIMUS_RAPHANgfnff_GFNFFca_structure)
        SEToriginal_SETraphanca_diffs = atom_coordinates_diferences(original_structure, PROPTIMUS_RAPHANgfnff_GFNFFca_structure)

        SETraphan_SETca.extend(SETraphan_SETca_diffs)
        SEToriginal_SETca.extend(SEToriginal_SETca_diffs)
        SETca_SETraphanca.extend(SETca_SETraphanca_diffs)
        SETraphan_SETraphanca.extend(SETraphan_SETraphanca_diffs)
        SEToriginal_SETraphanca.extend(SEToriginal_SETraphanca_diffs)

        protein = biotite.load_structure(file_path=glob(f"{calculation_dir}/optimised_PDB/*optimised.pdb")[0],
                                         model=1)
        h_bonds = set(int(atom) for h_bond in struc.hbond(protein) for atom in h_bond[1:3])
        for atom, SETraphan_SETca_diff, SEToriginal_SETca_diff, SETca_SETraphanca_diff, SETraphan_SETraphanca_diff, SEToriginal_SETraphanca_diff in zip(PROPTIMUS_RAPHANgfnff_structure.get_atoms(),
                                                                                                                                                        SETraphan_SETca_diffs,
                                                                                                                                                        SEToriginal_SETca_diffs,
                                                                                                                                                        SETca_SETraphanca_diffs,
                                                                                                                                                        SETraphan_SETraphanca_diffs,
                                                                                                                                                        SEToriginal_SETraphanca_diffs):
            atom._SETraphan_SETca_diff = SETraphan_SETca_diff
            atom._SEToriginal_SETca_diff = SEToriginal_SETca_diff
            atom._SETca_SETraphanca_diff = SETca_SETraphanca_diff
            atom._SETraphan_SETraphanca_diff = SETraphan_SETraphanca_diff
            atom._SEToriginal_SETraphanca_diff = SEToriginal_SETraphanca_diff

            if atom.name in ["H", "C", "N", "O", "CA"]:
                atom._h_bond = False
            else:
                atom._h_bond = atom.get_serial_number() - 1 in h_bonds
        structures.append(PROPTIMUS_RAPHANgfnff_structure)

    except FileNotFoundError:  # some optimisations did not converge
        continue
    except IndexError:
        continue
    counter += 1


print("Comparison of atomic positions:\n")
print(f"{counter} structures processed")
print(f"{len(SEToriginal_SETraphanca)} atoms processed\n")
print("                                       Mean, 50th percentile, 75th percentile, 90th percentile, 99th percentile, maximum value")
print(f"Statistics(SETraphan,          SETca): {statistics(SETraphan_SETca)}")
print(f"Statistics(SEToriginal,        SETca): {statistics(SEToriginal_SETca)}")
print("--------------------------------------------------------------------------------")
print(f"Statistics(SETraphan,   SETraphan+ca): {statistics(SETraphan_SETraphanca)}")
print(f"Statistics(SETca,       SETraphan+ca): {statistics(SETca_SETraphanca)}")
print(f"Statistics(SEToriginal, SETraphan+ca): {statistics(SEToriginal_SETraphanca)}")

create_histogram(SETraphan_SETca, "Atomic position deviations between  SET<sub>RAPHAN</sub> and SET<sub>Cα</sub>", 'Deviation of atomic positions [Å]', 'Number of atoms', "SETraphan_SETca_atomic_positions.pdf")
create_histogram(SEToriginal_SETca, "Atomic position deviations between SET<sub>ORIG</sub> and SET<sub>Cα</sub>", 'Deviation of atomic positions [Å]', 'Number of atoms', "SEToriginal_SETca_atomic_positions.pdf")
create_histogram(SETca_SETraphanca, "Atomic position deviations between SET<sub>Cα</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of atomic positions [Å]', 'Number of atoms', "SETca_SETraphan+ca_atomic_positions.pdf")
create_histogram(SETraphan_SETraphanca, "Atomic position deviations between SET<sub>RAPHAN</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of atomic positions [Å]', 'Number of atoms', "SETraphan_SETraphan+ca_atomic_positions.pdf")
create_histogram(SEToriginal_SETraphanca, "Atomic position deviations between SET<sub>ORIG</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of atomic positions [Å]', 'Number of atoms', "SEToriginal_SETraphan+ca_atomic_positions.pdf")

create_histogram_hbonds(structures, "_SETraphan_SETca_diff", SETraphan_SETca, 'SETraphan_SETca_atomic_positions_Hbonds.pdf')
create_histogram_hbonds(structures, "_SEToriginal_SETca_diff", SEToriginal_SETca, 'SEToriginal_SETca_atomic_positions_Hbonds.pdf')
create_histogram_hbonds(structures, "_SETca_SETraphanca_diff", SETca_SETraphanca, 'SETca_SETraphan+ca_atomic_positions_Hbonds.pdf')
create_histogram_hbonds(structures, "_SETraphan_SETraphanca_diff", SETraphan_SETraphanca, 'SETraphan_SETraphan+ca_atomic_positions_Hbonds.pdf')
create_histogram_hbonds(structures, "_SEToriginal_SETraphanca_diff", SEToriginal_SETraphanca, 'SEToriginal_SETraphan+ca_atomic_positions_Hbonds.pdf')

create_histogram_hbonds(structures, "_SETraphan_SETca_diff", SETraphan_SETca, 'SETraphan_SETca_atomic_positions_Hbonds.pdf', 'Atomic position deviations between SET<sub>RAPHAN</sub> and  SET<sub>Cα</sub>')
create_histogram_hbonds(structures, "_SEToriginal_SETca_diff", SEToriginal_SETca, 'SEToriginal_SETca_atomic_positions_Hbonds.pdf', 'Atomic position deviations between SET<sub>ORIG</sub> and SET<sub>Cα</sub>')
create_histogram_hbonds(structures, "_SETca_SETraphanca_diff", SETca_SETraphanca, 'SETca_SETraphan+ca_atomic_positions_Hbonds.pdf', 'Atomic position deviations between SET<sub>Cα</sub> and SET<sub>RAPHAN+Cα</sub>')
create_histogram_hbonds(structures, "_SETraphan_SETraphanca_diff", SETraphan_SETraphanca, 'SETraphan_SETraphan+ca_atomic_positions_Hbonds.pdf', 'Atomic position deviations between SET<sub>RAPHAN</sub> and SET<sub>RAPHAN+Cα</sub>')
create_histogram_hbonds(structures, "_SEToriginal_SETraphanca_diff", SEToriginal_SETraphanca, 'SEToriginal_SETraphan+ca_atomic_positions_Hbonds.pdf', 'Atomic position deviations between SET<sub>ORIG</sub> and SET<sub>RAPHAN+Cα</sub>')
