import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
from medication_dictionary import medication_dictionary
import py3Dmol

st.set_page_config(page_title="AI Healthcare Companion", layout="centered")
st.title("🩺 AI Healthcare Companion")

# Fetch chemical info from PubChem
def fetch_chemical_info(drug_name):
    compounds = pcp.get_compounds(drug_name, 'name')
    if compounds:
        c = compounds[0]
        info = {
            "iupac_name": c.iupac_name,
            "molecular_formula": c.molecular_formula,
            "molecular_weight": c.molecular_weight,
            "smiles": c.isomeric_smiles
        }
        return info
    return None

# Floating 2D molecule with atom labels + formula
def show_floating_2d_molecule(smiles, formula=None):
    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)  # generate 2D coordinates
    mol_block = rdmolfiles.MolToMolBlock(mol)
    
    view = py3Dmol.view(width=450, height=450)
    view.addModel(mol_block, 'mol')
    view.setStyle({'stick': {'colorscheme': 'cyanCarbon'}})
    view.setBackgroundColor('0x000000')
    view.rotate(30, 'x')
    view.rotate(15, 'y')
    view.zoomTo()
    
    # Atom labels
    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        view.addLabel(atom.GetSymbol(),
                      {'position': {'x': pos.x, 'y': pos.y, 'z': pos.z},
                       'backgroundColor': 'black',
                       'fontColor': 'cyan',
                       'fontSize': 14,
                       'showBackground': True})
    
    # Molecular formula label
    if formula:
        x_center = sum([conf.GetAtomPosition(a.GetIdx()).x for a in mol.GetAtoms()])/mol.GetNumAtoms()
        y_top = max([conf.GetAtomPosition(a.GetIdx()).y for a in mol.GetAtoms()]) + 1.5
        z_center = sum([conf.GetAtomPosition(a.GetIdx()).z for a in mol.GetAtoms()])/mol.GetNumAtoms()
        view.addLabel(formula,
                      {'position': {'x': x_center, 'y': y_top, 'z': z_center},
                       'backgroundColor': 'black',
                       'fontColor': 'lime',
                       'fontSize': 18,
                       'showBackground': True})

    html = view._make_html()
    st.components.v1.html(html, height=500, width=500)

# Streamlit GUI
symptom = st.text_input("Enter your symptom:")

if symptom:
    drugs = medication_dictionary.get(symptom.lower())
    
    if drugs:
        st.subheader(f"Possible medications for '{symptom}':")
        
        for drug in drugs:
            st.write(f"### {drug['name']}")
            st.write("**Dosage:**", drug['dosage'])
            
            # Warnings
            if "warnings" in drug and drug['warnings']:
                st.write("**Warnings:**")
                for w in drug['warnings']:
                    st.write(f"- {w}")
            
            # Chemical info
            info = fetch_chemical_info(drug['name'])
            if info:
                st.write("**IUPAC Name:**", info['iupac_name'])
                st.write("**Molecular Formula:**", info['molecular_formula'])
                st.write("**Molecular Weight:**", info['molecular_weight'])
                
                # Floating 2D molecule with annotations
                try:
                    show_floating_2d_molecule(info['smiles'], formula=info['molecular_formula'])
                except:
                    st.write("Could not render molecule.")
            else:
                st.write("Chemical info not found.")
    else:
        st.write("No medications found for this symptom.")
