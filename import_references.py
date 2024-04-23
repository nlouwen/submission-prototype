#!/usr/bin/env python3
# populate db reference table from existing data files

import json
import csv
import re
import sys
from pathlib import Path

from submission import create_app
from submission.extensions import db
from submission.models import Reference, Entry
from submission.utils.custom_errors import ReferenceNotFound


def main():
    app = create_app()
    with app.app_context():
        for filepath in Path("data").glob("*.json"):
            try:
                load_references(filepath)
            except Exception as e:
                print(f"Error loading references for {filepath}: {e}", file=sys.stderr)
        load_pending()


def load_references(filename):
    with open(filename, "r") as inf:
        if raw := inf.read():
            data = json.loads(raw)
        else:
            print(filename, "has no references", file=sys.stderr)
            return

        refs = set()
        for section in data:
            for field_name, field_value in data[section]:
                if field_name.endswith("references") and field_value:
                    refs.update(next(csv.reader([field_value], skipinitialspace=True)))
        valid_refs = []
        for ref in refs:
            if ";" in ref:
                ref, rest = ref.split(";", 1)
                valid_refs.append(rest)
            ref = re.sub(r"('|\")+", "", ref)

            if ref.lower().startswith("doi:"):
                ref = "doi:" + ref[4:]
            elif ref.startswith("https://doi.org/"):
                ref = "doi:" + ref[16:]
            elif ref.startswith("doi.org/"):
                ref = "doi:" + ref[8:]
            elif ref.startswith("10.") and "/" in ref:
                ref = "doi:" + ref
            elif ref.isdigit():
                ref = "pubmed:" + ref
            elif ref.startswith("PMID:"):
                ref = "pubmed:" + ref[5:]
            elif ref.find("doi/") != -1:
                idx = ref.find("10.", ref.find("doi/"))
                if idx != -1:
                    ref = "doi:" + ref[idx:]
            elif ref.startswith("https://patents.google.com/patent/"):
                ref = "patent:" + ref[33:]
            elif ref.lower().startswith("pubmed:"):
                ref = "pubmed:" + ref[7:]
            elif ref.startswith("https://"):
                ref = "url:" + ref

            ref = re.sub(r"\s+", "", ref)
            ref = re.sub(r"%2F", "/", ref)

            if valid_format(ref):
                valid_refs.append(ref)
            else:
                print(f"{filename}: Skipping invalid ref format: {ref!r}", file=sys.stderr)
        try:
            references = Reference.load_missing(valid_refs)
        except ReferenceNotFound as e:
            print(f"{filename}: {e}", file=sys.stderr)
            return

        bgc_id = filename.stem.replace("_data", "")
        entry = Entry.get_or_create(bgc_id=bgc_id)
        entry.add_references(references)


def valid_format(ref):
    regex = r"^pubmed:(\d+)$|^doi:10\.\d{4,9}/[-\\._;()/:a-zA-Z0-9]+$|^patent:(.+)$|^url:https?:\/\/(www\.)?[-a-zA-Z0-9@:%._\+~#=]{2,256}\.[a-z]{2,6}\b([-a-zA-Z0-9@:%_\+.~#?&//=]*)$"
    if re.match(regex, ref):
        return True
    return False


def load_pending():
    """Load the 'pending' keyword reference"""
    if not Reference.get("doi:pending"):
        ref = Reference(doi="pending")
        db.session.add(ref)
        db.session.commit()


if __name__ == "__main__":
    main()
