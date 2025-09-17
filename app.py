import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from medication_dictionary import medication_dictionary

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
                
                # Render 2D structure
                try:
                    mol = Chem.MolFromSmiles(info['smiles'])
                    img = Draw.MolToImage(mol)
                    st.image(img, caption=f"2D structure of {drug['name']}")
                except:
                    st.write("Could not render chemical structure.")
            else:
                st.write("Chemical info not found.")
    else:
        st.write("No medications found for this symptom.")
