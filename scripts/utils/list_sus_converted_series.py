#!/usr/bin/env python3
"""Enumerate converted SUS series per patient/timepoint and export to CSV."""
import csv
from pathlib import Path

BASES = {
    "BL": Path("/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_BL_converted_to_nifti"),
    "3Y": Path("/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_3Y_converted_to_nifti"),
    "5Y": Path("/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/data/SUS_conversion/SUS_5Y_converted_to_nifti"),
}

OUTPUT = Path("/nfs/br1_prosjekt/ParkWest/user/2026vae/AmundEspen/pd_thesis/scripts/sus_converted_series.csv")


def iter_series(base: Path, patient_dir: Path):
    if not patient_dir.exists():
        return []
    # list directories that contain any .nii* file
    series = []
    # include files directly under patient_dir
    direct_files = list(patient_dir.glob("*.nii*"))
    if direct_files:
        series.append((Path("."), direct_files))
    for subdir in sorted(p for p in patient_dir.iterdir() if p.is_dir()):
        nifti_files = list(subdir.glob("*.nii*"))
        if not nifti_files:
            for nested in sorted(p for p in subdir.iterdir() if p.is_dir()):
                nifti_files2 = list(nested.glob("*.nii*"))
                if nifti_files2:
                    series.append((nested.relative_to(patient_dir), nifti_files2))
            continue
        series.append((subdir.relative_to(patient_dir), nifti_files))
    return series


def main():
    rows = []
    for tp, base in BASES.items():
        if not base.exists():
            continue
        for patient_dir in sorted(p for p in base.iterdir() if p.is_dir()):
            patient = patient_dir.name
            series_list = iter_series(base, patient_dir)
            if not series_list:
                rows.append([tp, patient, "(none)", "0", ""])
                continue
            for relpath, files in series_list:
                file_names = "; ".join(f.name for f in files)
                rows.append([tp, patient, str(relpath), str(len(files)), file_names])

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["timepoint", "patient", "series_path", "n_files", "files"])
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {OUTPUT}")


if __name__ == "__main__":
    main()
