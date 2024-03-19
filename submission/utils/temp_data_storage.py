""" Temporary utils for data storage into a rough json file """

import json
from pathlib import Path
import subprocess
import datetime

from flask import current_app
from werkzeug.datastructures import MultiDict


class Storage:
    """Container class for storage methods"""

    @staticmethod
    def save_data(
        bgc_id: str, section_key: str, req_data: MultiDict, current_user
    ) -> None:
        """Append data to file containing all answers for one BGC

        Args:
            bgc_id (str): BGC identifier
            section_key (str): Section which form was filled in
            req_data (MultiDict): Submitted data coming from the request
        """
        data_dir = Path(current_app.root_path).parent / "data"

        data = {section_key: [(k, v) for k in req_data for v in req_data.getlist(k)]}
        existing_data = Storage.read_data(bgc_id)

        existing_data.update(data)

        if not existing_data.get("Changelog"):
            existing_data["Changelog"] = []
        existing_data["Changelog"].append(
            {
                "Edited_by": current_user.id,
                "Edited_at": datetime.datetime.now(datetime.UTC).strftime(
                    "%m/%d/%Y, %H:%M:%S"
                ),
            }
        )
        
        filename =  data_dir / f"{bgc_id}_data.json"
        with open(filename, "w") as outf:
            json.dump(existing_data, outf, sort_keys=True, indent=4)

        subprocess.run(f"git add {filename} && git commit -m 'updating {filename.name}'", shell=True, cwd=data_dir)

    @staticmethod
    def read_data(bgc_id: str) -> dict:
        """Read bgc data from an existing json file

        Args:
            bgc_id (str): BGC identifier

        Returns:
            dict: mapping of attributes to values
        """
        data_dir = Path(current_app.root_path).parent / "data"

        with open(data_dir / f"{bgc_id}_data.json", "r") as inf:
            if content := inf.read():
                existing_data = json.loads(content)
            else:
                existing_data = {}
        return existing_data

    @staticmethod
    def create_new_entry() -> str:
        """Look for existing files and create new entry

        Returns:
            str: new BGC identifier
        """
        data_dir = Path(current_app.root_path).parent / "data"

        max_entry_id = 0
        # TODO: create temp IDs
        for filepath in data_dir.glob("new*_data.json"):
            if (nr := int(filepath.stem[3:6])) > max_entry_id:
                max_entry_id = nr
        bgc_id = f"new{max_entry_id+1:0>3}"
        Path(data_dir / f"{bgc_id}_data.json").touch()

        return bgc_id

    @staticmethod
    def create_entry_if_not_exists(bgc_id: str):
        """Create data file for existing entry if it is somehow missing

        Args:
            bgc_id (str): MIBiG BGC identifier
        """
        data_dir = Path(current_app.root_path).parent / "data"

        if not (fname := data_dir / f"{bgc_id}_data.json").exists():
            fname.touch()
