""" Collection of custom validators used throughout the submission system """

import re
from wtforms import ValidationError
from pathlib import Path


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
    valid_ids = [f"BGC{num:0>7}" for num in range(1, 2750)]
    if bgc_id in valid_ids:
        return True
    if Path(f"{bgc_id}_data.json").exists():
        return True
    return False


# TODO: db cross validator
