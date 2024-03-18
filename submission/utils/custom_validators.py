""" Collection of custom validators used throughout the submission system """

import re
from pathlib import Path

from flask import current_app
from wtforms import ValidationError, validators


# TODO: rework to validate if entered gene IDs fall in locus
class ValidateTagListRegexp(object):
    def __init__(self, regex, message="Invalid Gene ID(s)"):
        self.regex = regex
        self.message = message

    def __call__(self, form, field):
        for gene_id in field.data:
            if not re.match(self.regex, gene_id):
                raise ValidationError(self.message)


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


# TODO: db cross validator
