from typing import Union
from pathlib import Path
from csv import DictReader

from sqlalchemy import select, inspect
from sqlalchemy.orm import Mapped, mapped_column

from submission.extensions import db


class NPAtlas(db.Model):
    """Static collection of NPAtlas entries"""

    __tablename__ = "npatlas"
    __table_args__ = {"schema": "edit"}

    id: Mapped[int] = mapped_column(autoincrement=True, primary_key=True)
    npaid: Mapped[str]
    compound_names: Mapped[str]
    compound_molecular_formula: Mapped[str]
    compound_accurate_mass: Mapped[float]
    compound_inchi: Mapped[str]
    compound_inchikey: Mapped[str]
    compound_smiles: Mapped[str]
    npatlas_url: Mapped[str]

    @staticmethod
    def get(compound: str) -> Union["NPAtlas", None]:
        """Get an NPAtlas entry by compound name

        Args:
            compound (str): compound name

        Returns:
            "NPAtlas" | None: NPAtlas entry or None if not exists
        """
        return db.session.scalar(
            select(NPAtlas).where(NPAtlas.compound_names == compound)
        )

    @staticmethod
    def from_tsv_dump(filepath: Path) -> None:
        """Fill the NPAtlas table from a tsv dump file

        Args:
            filepath (Path): path to the tsv dump file
        """
        properties = inspect(NPAtlas).attrs
        with open(filepath) as inf:
            data = DictReader(inf, dialect="excel-tab")
            for row in data:
                npa_entry = NPAtlas(**{k: v for k, v in row.items() if k in properties})
                db.session.add(npa_entry)
        db.session.commit()
