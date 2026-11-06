"""Download public VACV activity records from ChEMBL, preserving measurement types.

The activity table is a discovery pool, NOT extra SI training labels.
"""
import csv
import json
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

TARGETS = ("CHEMBL613493", "CHEMBL613168")  # Vaccinia virus; Vaccinia virus WR
BASE = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
FIELDS = (
    "activity_id", "target_chembl_id", "target_pref_name", "molecule_chembl_id",
    "canonical_smiles", "standard_type", "standard_relation", "standard_value",
    "standard_units", "assay_chembl_id", "assay_description", "document_chembl_id",
    "document_year", "data_validity_comment", "potential_duplicate", "source_url",
)


def fetch(target):
    offset = 0
    while True:
        query = urllib.parse.urlencode({"target_chembl_id": target, "limit": 1000, "offset": offset})
        with urllib.request.urlopen(f"{BASE}?{query}", timeout=60) as response:
            page = json.load(response)
        yield from page["activities"]
        offset += len(page["activities"])
        if not page["page_meta"]["next"]:
            break


def main():
    output = Path(__file__).resolve().parents[1] / "external"
    output.mkdir(exist_ok=True)
    records = {str(row["activity_id"]): row for target in TARGETS for row in fetch(target)}
    rows = sorted(records.values(), key=lambda row: int(row["activity_id"]))
    with (output / "chembl_vacv_activities.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=FIELDS)
        writer.writeheader()
        for record in rows:
            writer.writerow({key: (f"https://www.ebi.ac.uk/chembl/compound_report_card/{record['molecule_chembl_id']}/" if key == "source_url" else record.get(key)) for key in FIELDS})
    summary = {
        "source": "ChEMBL API", "targets": TARGETS, "total_activity_rows": len(rows),
        "unique_molecules": len({row["molecule_chembl_id"] for row in rows}),
        "measurement_types": dict(Counter(row["standard_type"] for row in rows)),
        "caution": "Activity rows are not interchangeable with SI labels. Do not add to SI training data without matching IC50 and CC50, compatible assay conditions, validated structures, and competition leakage checks.",
    }
    (output / "chembl_vacv_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Saved {len(rows)} activity rows from {len(TARGETS)} targets")


if __name__ == "__main__":
    main()
