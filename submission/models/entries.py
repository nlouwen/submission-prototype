from typing import Any, Union

from sqlalchemy import select
from sqlalchemy.orm import Mapped, mapped_column, relationship

from submission.extensions import db
from submission.models import Reference


# TODO: expand table
class Entry(db.Model):
    __tablename__ = "entries"
    __table_args__ = {"schema": "edit"}

    id: Mapped[int] = mapped_column(autoincrement=True, primary_key=True)
    identifier: Mapped[str]
    references: Mapped[list["Reference"]] = relationship(
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

    @staticmethod
    def get(bgc_id: str) -> Union["Entry", None]:
        """Get an entry from database based on identifier

        Args:
            bgc_id (str): BGC identifier

        Returns:
            Entry | None: entry database object or none if not exists
        """
        return db.session.scalar(select(Entry).where(Entry.identifier == bgc_id))

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
        entry = Entry.get(bgc_id=bgc_id)
        if entry is None:
            entry = Entry.create(bgc_id=bgc_id)

        refs = set()
        for locus in data["loci"]:
            for evidence in locus["evidence"]:
                refs.update(evidence["references"])

        loaded_refs = Reference.load_missing(list(refs))
        entry.references = loaded_refs

        db.session.commit()

    # TODO: create batch entries from data dir
