""" Collection of custom validators used throughout the submission system """

import re
from pathlib import Path

from flask import current_app
from rdkit import Chem
from wtforms import ValidationError, validators


# TODO: rework to validate if entered gene IDs fall in locus
class ValidateTagListRegexp(object):
    def __init__(self, regex, message="Invalid Gene ID: "):
        self.regex = regex
        self.message = message

    def __call__(self, form, field):
        for member in field.data:
            if not re.match(self.regex, member):
                message = self.message + member
                raise ValidationError(message)


class ValidateCitations(ValidateTagListRegexp):
    def __init__(self):
        regex = r"^doi:pending$|^pubmed:(\d+)$|^doi:10\.\d{4,9}/[-\\._;()/:a-zA-Z0-9]+$|^patent:(.+)$|^url:https?:\/\/(www\.)?[-a-zA-Z0-9@:%._\+~#=]{2,256}\.[a-z]{2,6}\b([-a-zA-Z0-9@:%_\+.~#?&//=]*)$"
        message = "Invalid citation format: "
        super().__init__(regex, message)


class ValidateStructureCross(ValidateTagListRegexp):
    def __init__(self):
        regex = r"^pubchem:(\d+)$|^chebi:(\d+)$|^chembl:CHEMBL(\d+)$|^chemspider:(\d+)$|^npatlas:NPA(\d+)$|^lotus:Q(\d+)$|^gnps:MSV(\d+)$|^cyanometdb:CyanoMetDB_(\d{4,4})$"
        message = "Invalid cross-reference format: "
        super().__init__(regex, message)


class ValidateEnzymeCross(ValidateTagListRegexp):
    def __init__(self):
        regex = r"^uniprot:[A-Z0-9]+$|^genpept:[A-Z]{3}[0-9]{5,7}\.[0-9]+$"
        message = "Invalid cross-reference format: "
        super().__init__(regex, message)


class ValidateReactionCross(ValidateTagListRegexp):
    def __init__(self):
        regex = r"^rhea:(\d+)$|^MITE(\d{7,7})$|^EC [0-9]+\.[0-9]+\.[0-9]+\.[0-9]+$"
        message = "Invalid cross-reference format: "
        super().__init__(regex, message)


def is_valid_bgc_id(bgc_id: str):
    valid_ids = [f"BGC{num:0>7}" for num in range(1, 3000)]
    if bgc_id in valid_ids:
        return True

    data_dir = Path(current_app.root_path).parent / "data"
    if Path(data_dir / f"{bgc_id}_data.json").exists():
        return True
    return False


class RequiredIf(validators.InputRequired):
    """Input required validator only if another field is filled

    Arguments:
        other_field_name (str): name of field to check
    """

    def __init__(self, other_field_name, *args, **kwargs):
        self.other_field_name = other_field_name
        super(RequiredIf, self).__init__(*args, **kwargs)
        self.field_flags.pop("required")

    def __call__(self, form, field):
        other_field = form._fields.get(self.other_field_name)
        if other_field is None:
            raise Exception(f"no field named {self.other_field_name} in form")
        if bool(other_field.data):
            super(RequiredIf, self).__call__(form, field)


def validate_smiles(form, field):
    mol = Chem.MolFromSmiles(field.data)
    if mol is None:
        raise ValidationError("Invalid SMILES")
