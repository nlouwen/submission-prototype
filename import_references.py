#!/usr/bin/env python3
# populate db reference table from existing data files

import json
import csv
import re
from pathlib import Path

from submission import create_app
from submission.extensions import db
from submission.models import Reference, Entry
from submission.utils.custom_errors import ReferenceNotFound


def main():
    app = create_app()
    with app.app_context():
        for filepath in Path("data").glob("*.json"):
            load_references(filepath)
        load_pending()


def load_references(filename):
    with open(filename, "r") as inf:
        if raw := inf.read():
            data = json.loads(raw)
        else:
            print(filename, "has no references")
            return

        refs = set()
        for section in data:
            for field_name, field_value in data[section]:
                if field_name.endswith("references") and field_value:
                    refs.update(next(csv.reader([field_value], skipinitialspace=True)))
        valid_refs = []
        for ref in refs:
            if valid_format(ref):
                valid_refs.append(ref)
            else:
                print(f"{filename}: Skipping invalid ref format: {ref}")
        try:
            references = Reference.load_missing(valid_refs)
        except ReferenceNotFound as e:
            raise Exception(str(filename) + ": " + str(e))

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
