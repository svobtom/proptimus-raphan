from glob import glob
from pathlib import Path

import plotly.graph_objects as go
import plotly.io as pio
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.dihedrals import Dihedral


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


def compare_dihedrals(file1, file2, d_phis, d_psis, d_omegas, d_chi1s):
    protein1 = mda.Universe(file1)
    protein2 = mda.Universe(file2)
    protein1_phis = Dihedral([res.phi_selection() for res in protein1.residues[1:]]).run().angles[0]
    protein1_psis = Dihedral([res.psi_selection() for res in protein1.residues[:-1]]).run().angles[0]
    protein1_omegas = Dihedral([res.omega_selection() for res in protein1.residues[:-1]]).run().angles[0]
    protein1_chi1s = Dihedral([res.chi1_selection() for res in protein1.residues if res.chi1_selection()]).run().angles[0]
    protein2_phis = Dihedral([res.phi_selection() for res in protein2.residues[1:]]).run().angles[0]
    protein2_psis = Dihedral([res.psi_selection() for res in protein2.residues[:-1]]).run().angles[0]
    protein2_omegas = Dihedral([res.omega_selection() for res in protein2.residues[:-1]]).run().angles[0]
    protein2_chi1s = Dihedral([res.chi1_selection() for res in protein2.residues if res.chi1_selection()]).run().angles[0]
    phis_diffs = np.abs(protein1_phis-protein2_phis)
    psis_diffs = np.abs(protein1_psis-protein2_psis)
    omegas_diffs = np.abs(protein1_omegas-protein2_omegas)
    chi1s_diffs = np.abs(protein1_chi1s-protein2_chi1s)
    phis_diffs = np.where(phis_diffs>180, 360-phis_diffs, phis_diffs)
    psis_diffs = np.where(psis_diffs>180, 360-psis_diffs, psis_diffs)
    omegas_diffs = np.where(omegas_diffs>180, 360-omegas_diffs, omegas_diffs)
    chi1s_diffs = np.where(chi1s_diffs>180, 360-chi1s_diffs, chi1s_diffs)
    d_phis.append(phis_diffs)
    d_psis.append(psis_diffs)
    d_omegas.append(omegas_diffs)
    d_chi1s.append(chi1s_diffs)




d_phis = []
d_psis = []
d_omegas = []
d_chi1s = []
counter = 0
for dir in sorted(glob("../calculations/AF-*")):
    if not Path(f"{dir}/constrained_alpha_carbons_optimisations/original/xtbopt.pdb").exists():
        continue
    compare_dihedrals(glob(f"{dir}/constrained_alpha_carbons_optimisations/original/xtbopt.pdb")[0], glob(f"{dir}/optimised_PDB/*optimised.pdb")[0], d_phis, d_psis, d_omegas, d_chi1s)
    counter += 1
d_phis = np.concatenate(d_phis)
d_psis = np.concatenate(d_psis)
d_omegas = np.concatenate(d_omegas)
d_chi1s = np.concatenate(d_chi1s)
print("\n\nComparison of selected dihedral angles between SETraphan and SETca")
print(f"{counter} structures processed")
print("                     Mean, 50th percentile, 75th percentile, 90th percentile, 99th percentile, maximum value")
print(f"Statistics(Phi)  : {statistics(d_phis)}")
print(f"Statistics(Psi)  : {statistics(d_psis)}")
print(f"Statistics(Omega): {statistics(d_omegas)}")
print(f"Statistics(Chi1) : {statistics(d_chi1s)}")
create_histogram(d_phis, "Phi dihedral angle deviations between SET<sub>RAPHAN</sub> and SET<sub>Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SETraphan_SETca_phi_dihedral_angles.pdf")
create_histogram(d_psis, "Psi dihedral angle deviations between SET<sub>RAPHAN</sub> and SET<sub>Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SETraphan_SETca_psi_dihedral_angles.pdf")
create_histogram(d_omegas, "Omega dihedral angle deviations between SET<sub>RAPHAN</sub> and SET<sub>Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SETraphan_SETca_omega_dihedral_angles.pdf")
create_histogram(d_chi1s, "Chi1 dihedral angle deviations between SET<sub>RAPHAN</sub> and SET<sub>Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SETraphan_SETca_chi1_dihedral_angles.pdf")



d_phis = []
d_psis = []
d_omegas = []
d_chi1s = []
counter = 0
for dir in sorted(glob("../calculations/AF-*")):
    if not Path(f"{dir}/constrained_alpha_carbons_optimisations/original/xtbopt.pdb").exists():
        continue
    compare_dihedrals(glob(f"{dir}/input_PDB/*")[0], f"{dir}/constrained_alpha_carbons_optimisations/original/xtbopt.pdb", d_phis, d_psis, d_omegas, d_chi1s)
    counter += 1
d_phis = np.concatenate(d_phis)
d_psis = np.concatenate(d_psis)
d_omegas = np.concatenate(d_omegas)
d_chi1s = np.concatenate(d_chi1s)
print("\nComparison of selected dihedral angles between SEToriginal and SETca")
print(f"{counter} structures processed")
print("                     Mean, 50th percentile, 75th percentile, 90th percentile, 99th percentile, maximum value")
print(f"Statistics(Phi)  : {statistics(d_phis)}")
print(f"Statistics(Psi)  : {statistics(d_psis)}")
print(f"Statistics(Omega): {statistics(d_omegas)}")
print(f"Statistics(Chi1) : {statistics(d_chi1s)}")
create_histogram(d_phis, "Phi dihedral angle deviations between SET<sub>ORIG</sub> and SET<sub>Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SEToriginal_SETca_phi_dihedral_angles.pdf")
create_histogram(d_psis, "Psi dihedral angle deviations between SET<sub>ORIG</sub> and SET<sub>Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SEToriginal_SETca_psi_dihedral_angles.pdf")
create_histogram(d_omegas, "Omega dihedral angle deviations between SET<sub>ORIG</sub> and SET<sub>Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SEToriginal_SETca_omega_dihedral_angles.pdf")
create_histogram(d_chi1s, "Chi1 dihedral angle deviations between SET<sub>ORIG</sub> and SET<sub>Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SEToriginal_SETca_chi1_dihedral_angles.pdf")



d_phis = []
d_psis = []
d_omegas = []
d_chi1s = []
counter = 0
for dir in sorted(glob("../calculations/AF-*")):
    if not Path(f"{dir}/constrained_alpha_carbons_optimisations/original/xtbopt.pdb").exists():
        continue
    compare_dihedrals(glob(f"{dir}/optimised_PDB/*optimised.pdb")[0], f"{dir}/constrained_alpha_carbons_optimisations/raphan/xtbopt.pdb", d_phis, d_psis, d_omegas, d_chi1s)
    counter += 1
d_phis = np.concatenate(d_phis)
d_psis = np.concatenate(d_psis)
d_omegas = np.concatenate(d_omegas)
d_chi1s = np.concatenate(d_chi1s)
print("\n\nComparison of selected dihedral angles between SETraphan and SETraphan+ca")
print(f"{counter} structures processed")
print("                     Mean, 50th percentile, 75th percentile, 90th percentile, 99th percentile, maximum value")
print(f"Statistics(Phi)  : {statistics(d_phis)}")
print(f"Statistics(Psi)  : {statistics(d_psis)}")
print(f"Statistics(Omega): {statistics(d_omegas)}")
print(f"Statistics(Chi1) : {statistics(d_chi1s)}")
create_histogram(d_phis, "Phi dihedral angle deviations between SET<sub>RAPHAN</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SETraphan_SETraphan+ca_phi_dihedral_angles.pdf")
create_histogram(d_psis, "Psi dihedral angle deviations between SET<sub>RAPHAN</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SETraphan_SETraphan+ca_psi_dihedral_angles.pdf")
create_histogram(d_omegas, "Omega dihedral angle deviations between SET<sub>RAPHAN</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SETraphan_SETraphan+ca_omega_dihedral_angles.pdf")
create_histogram(d_chi1s, "Chi1 dihedral angle deviations between SET<sub>RAPHAN</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SETraphan_SETraphan+ca_chi1_dihedral_angles.pdf")



d_phis = []
d_psis = []
d_omegas = []
d_chi1s = []
counter = 0
for dir in sorted(glob("../calculations/AF-*")):
    if not Path(f"{dir}/constrained_alpha_carbons_optimisations/original/xtbopt.pdb").exists():
        continue
    compare_dihedrals(glob(f"{dir}/constrained_alpha_carbons_optimisations/original/xtbopt.pdb")[0], f"{dir}/constrained_alpha_carbons_optimisations/raphan/xtbopt.pdb", d_phis, d_psis, d_omegas, d_chi1s)
    counter += 1
d_phis = np.concatenate(d_phis)
d_psis = np.concatenate(d_psis)
d_omegas = np.concatenate(d_omegas)
d_chi1s = np.concatenate(d_chi1s)
print("\n\nComparison of selected dihedral angles between SETca and SETraphanca")
print(f"{counter} structures processed")
print("                     Mean, 50th percentile, 75th percentile, 90th percentile, 99th percentile, maximum value")
print(f"Statistics(Phi)  : {statistics(d_phis)}")
print(f"Statistics(Psi)  : {statistics(d_psis)}")
print(f"Statistics(Omega): {statistics(d_omegas)}")
print(f"Statistics(Chi1) : {statistics(d_chi1s)}")
create_histogram(d_phis, "Phi dihedral angle deviations between SET<sub>Cα</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SETca_SETraphan+ca_phi_dihedral_angles.pdf")
create_histogram(d_psis, "Psi dihedral angle deviations between SET<sub>Cα</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SETca_SETraphan+ca_psi_dihedral_angles.pdf")
create_histogram(d_omegas, "Omega dihedral angle deviations between SET<sub>Cα</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SETca_SETraphan+ca_omega_dihedral_angles.pdf")
create_histogram(d_chi1s, "Chi1 dihedral angle deviations between SET<sub>Cα</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SETca_SETraphan+ca_chi1_dihedral_angles.pdf")



d_phis = []
d_psis = []
d_omegas = []
d_chi1s = []
counter = 0
for dir in sorted(glob("../calculations/AF-*")):
    if not Path(f"{dir}/constrained_alpha_carbons_optimisations/original/xtbopt.pdb").exists():
        continue
    compare_dihedrals(glob(f"{dir}/input_PDB/*")[0], f"{dir}/constrained_alpha_carbons_optimisations/raphan/xtbopt.pdb", d_phis, d_psis, d_omegas, d_chi1s)
    counter += 1
d_phis = np.concatenate(d_phis)
d_psis = np.concatenate(d_psis)
d_omegas = np.concatenate(d_omegas)
d_chi1s = np.concatenate(d_chi1s)
print("\nComparison of selected dihedral angles between SEToriginal and SETraphanca")
print(f"{counter} structures processed")
print("                     Mean, 50th percentile, 75th percentile, 90th percentile, 99th percentile, maximum value")
print(f"Statistics(Phi)  : {statistics(d_phis)}")
print(f"Statistics(Psi)  : {statistics(d_psis)}")
print(f"Statistics(Omega): {statistics(d_omegas)}")
print(f"Statistics(Chi1) : {statistics(d_chi1s)}")
create_histogram(d_phis, "Phi dihedral angle deviations between SET<sub>ORIG</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SEToriginal_SETraphan+ca_phi_dihedral_angles.pdf")
create_histogram(d_psis, "Psi dihedral angle deviations between SET<sub>ORIG</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SEToriginal_SETraphan+ca_psi_dihedral_angles.pdf")
create_histogram(d_omegas, "Omega dihedral angle deviations between SET<sub>ORIG</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SEToriginal_SETraphan+ca_omega_dihedral_angles.pdf")
create_histogram(d_chi1s, "Chi1 dihedral angle deviations between SET<sub>ORIG</sub> and SET<sub>RAPHAN+Cα</sub>", 'Deviation of dihedral angles [°]', 'Number of dihedral angles', "SEToriginal_SETraphan+ca_chi1_dihedral_angles.pdf")

