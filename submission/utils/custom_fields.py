""" Collection of custom form Field classes used throughout the submission system"""

import re
import csv
from wtforms import Field, StringField, SelectMultipleField, validators, widgets
from typing import Any, Optional

from submission.utils.custom_validators import validate_smiles
from submission.utils.custom_widgets import StructureInput


class GeneIdField(StringField):
    """Reusable Gene ID field with validators"""

    def __init__(
        self,
        label: Optional[str] = None,
        validators: Optional[list[Any]] = [
            validators.Optional(),
            validators.Regexp(r"^[^, ]*$", message="Invalid Gene ID"),
        ],
        description: Optional[
            str
        ] = "NCBI GenPept ID (e.g. CAB60185.1), locus tag (e.g. SCO6266), or gene name (e.g. scbA)",
        **kwargs,
    ):
        super(GeneIdField, self).__init__(
            label=label, validators=validators, description=description, **kwargs
        )


class TagListField(Field):
    """Custom field for comma separated input, use double quotes to enter names containing commas"""

    widget = widgets.TextInput()

    def _value(self):
        """convert list of tags into comma separated string"""
        if self.data:
            return '"%s"' % '", "'.join(self.data)
        else:
            return ""

    def process_formdata(self, valuelist: list[str]):
        """process comma separated string into list, remove duplicates

        Args:
            valuelist (list): form field input as ["..."]
        """
        if valuelist:
            unique_tags = set(next(csv.reader(valuelist, skipinitialspace=True)))
            self.data = list(unique_tags)
        else:
            self.data = []


class ReferenceField(TagListField):

    def process_formdata(self, valuelist: list[str]):
        """Process references, convert to standardized format, removes duplicates

        Args:
            valuelist (list): form field input as ["..."]
        """
        super().process_formdata(valuelist)

        sanitized_refs = set()
        for ref in self.data:
            ref = re.sub(r"('|\")+", "", ref)  # remove quotes
            ref = re.sub(r"\s+", "", ref)  # remove whitespaces
            ref = re.sub(r"%2F", "/", ref)  # revert html url encoding

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

            sanitized_refs.add(ref)
        self.data = list(sanitized_refs)


class MultiCheckboxField(SelectMultipleField):
    """
    A multiple-select, except displays a list of checkboxes.

    Iterating the field will produce subfields, allowing custom rendering of
    the enclosed checkbox fields.
    """

    widget = widgets.ListWidget(prefix_label=False)
    option_widget = widgets.CheckboxInput()


def smiles_field_factory(
    label: Optional[str] = None,
    description: Optional[str] = None,
    required: bool = False,
    show_structure: bool = True,
):
    """Creates a customized SMILES input filed

    Args:
        label (Optional[str]): field label
        description (Optional[str]): field description
        required (bool): flag to add the InputRequired validator
        show_structure (bool): flag to add the StructureInput widget

    Returns:
        SmilesField: customized SMILES input field
    """
    if description is None:
        description = "SMILES representation of the structure, preferentially isomeric"

    default_validators = [
        validate_smiles,
    ]
    if required:
        default_validators.append(validators.InputRequired())
    else:
        default_validators.append(validators.Optional())

    class SmilesField(StringField):
        """Standardized SMILES input field"""

        if show_structure:
            widget = StructureInput()

        def __init__(
            self,
            label: Optional[str] = label,
            validators: list[Any] = default_validators,
            description: Optional[str] = description,
            **kwargs,
        ):
            super(SmilesField, self).__init__(
                label=label, validators=validators, description=description, **kwargs
            )

    return SmilesField()
