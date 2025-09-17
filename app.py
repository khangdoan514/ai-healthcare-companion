import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw

from medication_dictionary import medication_dictionary

st.title("AI Healthcare Companion")

# User input
symptom = st.text_input("Enter your symptom:")

if symptom:
    drugs = medication_dictionary.get(symptom.lower())
    if drugs:
        st.write(f"Possible medications for '{symptom}':")
        for drug in drugs:
            st.write(f"- {drug['name']}: {drug['dosage']}")

            # Fetch chemical info
            compounds = pcp.get_compounds(drug['name'], 'name')
            if compounds:
                c = compounds[0]
                st.write("IUPAC Name:", c.iupac_name)
                st.write("Molecular Formula:", c.molecular_formula)
                st.write("Molecular Weight:", c.molecular_weight)

                # Render 2D structure
                mol = Chem.MolFromSmiles(c.isomeric_smiles)
                img = Draw.MolToImage(mol)
                st.image(img, caption=f"2D structure of {drug['name']}")
            else:
                st.write("No chemical data found.")
    else:
        st.write("No medication found for this symptom.")
