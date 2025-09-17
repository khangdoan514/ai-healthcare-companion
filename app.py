import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw

st.set_page_config(page_title="AI Healthcare Companion", layout="centered")
st.title("🩺 AI Healthcare Companion")

# -------------------------
# Step 1: Symptom → Drug Mapping
# -------------------------
def suggest_drugs(symptom):
    """Suggest possible medications for a given symptom."""
    mapping = {
        "headache": ["ibuprofen", "acetaminophen"],
        "fever": ["acetaminophen"],
        "cough": ["dextromethorphan"]
    }
    return mapping.get(symptom.lower(), [])

# -------------------------
# Step 2: Fetch Chemical Info
# -------------------------
def fetch_chemical_info(drug_name):
    """Fetch chemical information from PubChem."""
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
# Step 3: Streamlit GUI
# -------------------------
symptom = st.text_input("Enter your symptom:")

if symptom:
    suggested_drugs = suggest_drugs(symptom)
    
    if suggested_drugs:
        st.subheader(f"Possible medications for '{symptom}':")
        for drug in suggested_drugs:
            st.write(f"### {drug}")
            
            info = fetch_chemical_info(drug)
            if info:
                st.write("**IUPAC Name:**", info['iupac_name'])
                st.write("**Molecular Formula:**", info['molecular_formula'])
                st.write("**Molecular Weight:**", info['molecular_weight'])

                # 2D Structure Rendering
                try:
                    mol = Chem.MolFromSmiles(info['smiles'])
                    img = Draw.MolToImage(mol)
                    st.image(img, caption=f"2D structure of {drug}")
                except:
                    st.write("Could not render chemical structure.")
            else:
                st.write("Chemical info not found.")
    else:
        st.write("No medications found for this symptom.")
