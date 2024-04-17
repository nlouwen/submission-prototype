""" Class to handle reference standardization """

import csv
import requests

from .temp_data_storage import Storage


class ReferenceUtils:

    @staticmethod
    def collect_references(bgc_id: str) -> set[str]:
        """Collect all known references for this entry

        Args:
            bgc_id (str): BGC identifier

        Returns:
            set[str]: references
        """
        data = Storage.read_data(bgc_id)

        refs = set()
        for section in data:
            for field_name, field_value in data[section]:
                if field_name.endswith("references") and field_value:
                    refs.update(next(csv.reader([field_value], skipinitialspace=True)))
        return refs

    @staticmethod
    def get_reference_metadata(reference: str) -> dict[str, str]:
        """Collect metadata on reference from doi/PMID using liningtonlab's api

        Args:
            reference (str): reference formatted as 'doi:10....' | 'pubmed:73....'

        Returns:
            dict[str, str]: mapping of reference metadata
        """
        id_type, identifier = reference.split(":", 1)
        if id_type == "pubmed":
            id_type = "pmid"
        r = requests.post(
            f"https://litapi.liningtonlab.org/article/?{id_type}={identifier}"
        )
        return r.json()

    @staticmethod
    def append_ref(current: str, new_ref: str) -> str:
        """Append a reference to an existing input

        Args:
            current_refs (str): already filled references, double quoted comma separated
            new_ref (str): reference to be added

        Returns:
            str: new input value with appended reference
        """
        if not current:
            current = f'"{new_ref}"'
        else:
            current_refs = set(next(csv.reader([current], skipinitialspace=True)))
            if new_ref not in current_refs:
                current += f', "{new_ref}"'
        return current
