from typing import Any

from sqlalchemy.orm import Mapped, mapped_column

from submission.extensions import db
from submission.models import Reference


# TODO: expand table
class Entry(db.Model):
    __tablename__ = "entries"
    __table_args__ = {"schema": "edit"}

    id: Mapped[int] = mapped_column(autoincrement=True, primary_key=True)
    identifier: Mapped[str]
    references = db.relationship(
        "Reference",
        secondary="edit.entry_references",
        back_populates="entries",
        lazy="selectin",
    )

    @classmethod
    def create(cls, bgc_id: str):
        """Create a database entry for this BGC

        Args:
            bgc_id (str): BGC identifier
        """
        entry = cls(identifier=bgc_id)
        db.session.add(entry)
        db.session.commit()
        return entry

    # TODO: save all important data
    @staticmethod
    def save_minimal(bgc_id: str, data: dict[str, Any]):
        """Save minimal entry information

        Args:
            bgc_id (str): BGC identfier
            data (dict): Minimal information to save
        """
        # if minimal is the first section submitted, no entry will exist in db yet
        # TODO: create entries separately and consistently over all sections
        entry = Entry.query.filter(Entry.identifier == bgc_id).first()
        if not entry:
            entry = Entry.create(bgc_id=bgc_id)

        refs = []
        for locus in data["loci"]:
            for evidence in locus["evidence"]:
                refs.extend(evidence["references"])

        loaded_refs = Reference.load_missing(refs)
        entry.references = loaded_refs

        db.session.commit()

    # TODO: create batch entries from data dir
