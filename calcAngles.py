"""
pdb/txt file backbone atom torsion angles calculator
Bioinformatics application
Author: Khushdeep Brar
December 2024
LLM assistance utilized to debug
GPT Version 4.0
"""

import numpy as np

def parse_pdb(file_path):
    """
    Parse PDB file to extract atomic coordinates of backbone atoms.
    """
    residues = {}
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                resid = int(line[22:26].strip())
                chain = line[21].strip()
                x, y, z = map(float, (line[30:38], line[38:46], line[46:54]))
                if resid not in residues:
                    residues[resid] = {}
                residues[resid][atom_name] = (x, y, z)
    return residues

def vector(p1, p2):
    """Calculate vector between two points."""
    return np.array(p2) - np.array(p1)

def normalize(v):
    """Normalize a vector."""
    return v / np.linalg.norm(v)

def calculate_dihedral(p1, p2, p3, p4):
    """
    Calculate dihedral angle between four points in space.
    """
    b1 = -vector(p2, p1)
    b2 = vector(p2, p3)
    b3 = vector(p3, p4)

    n1 = np.cross(normalize(b1), normalize(b2))
    n2 = np.cross(normalize(b2), normalize(b3))
    m1 = np.cross(n1, normalize(b2))

    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return np.degrees(np.arctan2(y, x))

def calculate_phi_psi(residues, start, end):
    """
    Calculate phi and psi torsion angles for a range of residues.
    """
    torsion_angles = []
    residue_ids = sorted(residues.keys())
    for i, resid in enumerate(residue_ids):
        if start <= resid <= end:
            try:
                if resid - 1 in residues and resid + 1 in residues:
                    prev_res = residues[resid - 1]
                    current_res = residues[resid]
                    next_res = residues[resid + 1]
                    phi = calculate_dihedral(prev_res["C"], current_res["N"], current_res["CA"], current_res["C"])
                    psi = calculate_dihedral(current_res["N"], current_res["CA"], current_res["C"], next_res["N"])
                    torsion_angles.append((resid, phi, psi))
                else:
                    torsion_angles.append((resid, None, None))  # Missing neighbors
            except KeyError:
                torsion_angles.append((resid, None, None))  # Missing backbone atoms
    return torsion_angles

def main():
    pdb_file = "C:/Users/Khushdeep.SKYBEAST/Downloads/model.txt"  # Path
    start_res = 80
    end_res = 90

    # Parse the PDB file
    residues = parse_pdb(pdb_file)

    # Calculate torsion angles
    torsion_angles = calculate_phi_psi(residues, start_res, end_res)

    print("Residue   Phi       Psi")
    print("-----------------------")
    for resid, phi, psi in torsion_angles:
        phi_val = f"{phi:.2f}" if phi is not None else "N/A"
        psi_val = f"{psi:.2f}" if psi is not None else "N/A"
        print(f"{resid:<8}{phi_val:<10}{psi_val:<10}")

if __name__ == "__main__":
    main()
