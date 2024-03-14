""" Collection of custom form Field classes used throughout the submission system"""

import csv
from wtforms import Field, StringField, SelectMultipleField, validators, widgets
from typing import Any


class GeneIdField(StringField):
    """Reusable Gene ID field with validators"""

    def __init__(
        self,
        label: str | None = None,
        validators: list[Any] | None = [
            validators.Optional(),
            validators.Regexp(r"^[^, ]*$", message="Invalid Gene ID"),
        ],
        **kwargs,
    ):
        super(GeneIdField, self).__init__(label, validators, **kwargs)


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
