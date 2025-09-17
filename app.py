import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
import io
import base64
import json

st.set_page_config(page_title="AI Healthcare Companion", layout="centered")
st.title("AI Healthcare Companion 🩺")

# ---------- Symptom → Medicine Mapping ----------
# You can expand this JSON later or move to SQLite
symptom_to_meds = {
    "headache": {"name": "ibuprofen", "dosage": "200-400 mg every 4-6 hours"},
    "fever": {"name": "acetaminophen", "dosage": "500-1000 mg every 4-6 hours"},
    "cough": {"name": "dextromethorphan", "dosage": "10-20 mg every 4 hours"}
}

# ---------- User Input ----------
symptom = st.text_input("Enter your symptom:")

if symptom:
    med_info = symptom_to_meds.get(symptom.lower())
    if med_info:
        med_name = med_info["name"]
        st.subheader(f"Suggested Medicine: {med_name}")
        st.write(f"Recommended dosage: {med_info['dosage']}")

        # ---------- Fetch Chemical Info from PubChem ----------
        compounds = pcp.get_compounds(med_name, 'name')
        if compounds:
            c = compounds[0]
            st.write("**IUPAC Name:**", c.iupac_name)
            st.write("**Molecular Formula:**", c.molecular_formula)
            st.write("**Molecular Weight:**", round(c.molecular_weight, 2))

            # ---------- Render 2D Structure Using RDKit ----------
            mol = Chem.MolFromSmiles(c.isomeric_smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(300, 300))
                st.image(img, caption=f"2D structure of {med_name}")

        else:
            st.write("No chemical data found in PubChem.")
    else:
        st.warning("No suggestion available for this symptom. Try another one!")