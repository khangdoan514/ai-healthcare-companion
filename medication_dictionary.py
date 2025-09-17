# medication_dictionary.py

# Maps symptoms to possible medications with dosage and warnings
medication_dictionary = {
    "headache": [
        {
            "name": "ibuprofen",
            "dosage": "200-400 mg every 4-6h",
            "warnings": ["Avoid if allergic to NSAIDs", "Take with food to prevent stomach upset"]
        },
        {
            "name": "acetaminophen",
            "dosage": "500-1000 mg every 4-6h",
            "warnings": ["Do not exceed 4000 mg per day", "Avoid alcohol"]
        }
    ],
    "fever": [
        {
            "name": "acetaminophen",
            "dosage": "500-1000 mg every 4-6h",
            "warnings": ["Do not exceed 4000 mg per day", "Avoid alcohol"]
        }
    ],
    "cough": [
        {
            "name": "dextromethorphan",
            "dosage": "10-20 mg every 4h",
            "warnings": ["Do not exceed 120 mg per day", "Avoid in children under 6"]
        }
    ]
}
