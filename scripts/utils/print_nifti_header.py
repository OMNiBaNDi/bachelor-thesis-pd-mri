#!/usr/bin/env python3
"""Print the NIfTI header (and any embedded JSON extension) for a given file."""

import argparse
import json
import sys
from pathlib import Path

try:
    import nibabel as nib  # type: ignore
except ModuleNotFoundError:
    print("nibabel is required (pip install --user nibabel)", file=sys.stderr)
    sys.exit(1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("nifti", help="Path to the .nii or .nii.gz file to inspect")
    args = parser.parse_args()

    nifti_path = Path(args.nifti).expanduser()
    if not nifti_path.exists():
        print(f"File not found: {nifti_path}", file=sys.stderr)
        sys.exit(1)

    img = nib.load(str(nifti_path))
    header = img.header
    print(f"=== {nifti_path} ===")
    print("-- Core header --")
    print(header)

    descrip = header.get("descrip")
    if descrip is not None:
        if hasattr(descrip, "tobytes"):
            descrip = descrip.tobytes()
        if isinstance(descrip, bytes):
            descrip = descrip.decode("utf-8", errors="ignore")
        print(f"descrip: {descrip!r}")

    if header.extensions:
        print("-- Extensions --")
        for idx, ext in enumerate(header.extensions):
            payload = ext.get_content()
            if isinstance(payload, bytes):
                text = payload.decode("utf-8", errors="ignore").strip()
            else:
                text = str(payload).strip()
            print(f"Extension {idx} (ecode={ext.get_code()}):")
            if text.startswith("{"):
                try:
                    obj = json.loads(text)
                    print(json.dumps(obj, indent=2))
                except json.JSONDecodeError:
                    print(text)
            else:
                print(text)
    else:
        print("-- No header extensions present --")


if __name__ == "__main__":
    main()
