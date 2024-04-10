""" Collection of custom form Field classes used throughout the submission system"""

import csv
from wtforms import Field, StringField, SelectMultipleField, validators, widgets
from typing import Any, Optional

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
        if self.data:
            return '"%s"' % '", "'.join(self.data)
        else:
            return ""

    def process_formdata(self, valuelist):
        if valuelist:
            self.data = next(csv.reader(valuelist, skipinitialspace=True))
        else:
            self.data = []


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
        validators.Regexp(
            r"^[\[\]a-zA-Z0-9\@()=\/\\#+.%*-]+$", message="Invalid SMILES"
        )
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
