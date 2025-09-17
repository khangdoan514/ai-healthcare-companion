import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
from medication_dictionary import medication_dictionary
import py3Dmol

st.set_page_config(page_title="AI Healthcare Companion", layout="centered")
st.title("🩺 AI Healthcare Companion")

# -------------------------
# Fetch chemical info from PubChem
# -------------------------
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

# -------------------------
# Render floating 2D molecule (Big Hero 6 style)
# -------------------------
def show_floating_2d_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)  # generate 2D coordinates
    mol_block = rdmolfiles.MolToMolBlock(mol)
    
    view = py3Dmol.view(width=400, height=400)
    view.addModel(mol_block, 'mol')
    view.setStyle({'stick': {'colorscheme': 'cyanCarbon'}})  # neon color
    view.setBackgroundColor('0x000000')  # black background
    view.rotate(30, 'x')  # slight tilt for floating effect
    view.rotate(15, 'y')
    view.zoomTo()
    html = view._make_html()
    st.components.v1.html(html, height=450, width=450)

# -------------------------
# Streamlit GUI
# -------------------------
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
                
                # Floating 2D molecule
                try:
                    show_floating_2d_molecule(info['smiles'])
                except:
                    st.write("Could not render molecule.")
            else:
                st.write("Chemical info not found.")
    else:
        st.write("No medications found for this symptom.")
