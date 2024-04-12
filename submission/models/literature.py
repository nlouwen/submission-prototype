from typing import Union

from sqlalchemy import select, or_
from sqlalchemy.orm import Mapped, mapped_column

from submission.extensions import db
from submission.utils import ReferenceUtils


class Reference(db.Model):
    __tablename__ = "references"
    __table_args__ = {"schema": "edit"}

    id: Mapped[int] = mapped_column(autoincrement=True, primary_key=True)
    doi: Mapped[str] = mapped_column(nullable=True)
    pubmed: Mapped[str] = mapped_column(nullable=True)
    title: Mapped[str] = mapped_column(nullable=True)
    authors: Mapped[str] = mapped_column(nullable=True)
    year: Mapped[int] = mapped_column(nullable=True)
    journal: Mapped[str] = mapped_column(nullable=True)
    entries = db.relationship(
        "Entry",
        secondary="edit.entry_references",
        back_populates="references",
        lazy="selectin",
    )

    @property
    def identifier(self):
        """Access a reference identifier, prefers doi

        If somehow no identifier is present, returns None
        """
        if self.doi:
            return f"doi:{self.doi}"
        elif self.pubmed:
            return f"pubmed:{self.pubmed}"
        else:
            return None

    def short_authors(self):
        """Shorten the authors to the first et al."""
        if not self.authors:
            return ""

        all_authors = self.authors.split(";")
        first_author_names = all_authors[0].split(",")

        if len(first_author_names) == 2:
            surname, firstname = first_author_names
            names = f"{surname.strip()}, {firstname.strip()[0]}"
        else:
            names = f"{first_author_names[0].strip()}"

        if len(all_authors) > 1:
            others = " et al."
        else:
            others = ""
        return f"{names}{others}"

    def summarize(self):
        """Generate a one-line summary of the reference"""
        title = self.title if self.title else ""
        journal = self.journal if self.journal else ""
        year = self.year if self.year else ""
        identifier = self.identifier if self.identifier else ""
        return f"{title} {self.short_authors()} {journal}, {year}. {identifier}"

    @classmethod
    def load(cls, reference: str):
        """Loads a reference into the reference table

        Args:
            reference (str): reference formatted as doi:10... | pubmed:73...
        """
        metadata = ReferenceUtils.get_reference_metadata(reference)
        ref = cls(
            doi=metadata.get("doi"),
            pubmed=metadata.get("pmid"),
            title=metadata.get("title"),
            authors=metadata.get("authors"),
            year=metadata.get("year"),
            journal=metadata.get("journal"),
        )
        db.session.add(ref)
        db.session.commit()
        return ref

    @staticmethod
    def load_missing(references: list[str]) -> list["Reference"]:
        """Load missing references into database

        Args:
            references (list[str]): references to load

        Returns:
            list[Reference]: list of all references, either existing or newly loaded
        """
        refs = []
        for reference in references:
            ref = Reference.get(reference)
            if ref is None:
                ref = Reference.load(reference)
            refs.append(ref)
        return refs

    @staticmethod
    def get(reference: str) -> Union["Reference", None]:
        """Get a reference object from the database based on identifier

        Args:
            reference (str): reference formatted as doi:10... | pubmed:77...

        Returns:
            Reference | None: reference database object or None if not exists
        """
        ident = reference.split(":", 1)[1]
        return db.session.scalar(
            select(Reference).where(
                or_(Reference.doi == ident, Reference.pubmed == ident)
            )
        )


class EntryReference(db.Model):
    __tablename__ = "entry_references"
    __table_args__ = {"schema": "edit"}

    entry_id: Mapped[int] = mapped_column(
        db.ForeignKey("edit.entries.id"), primary_key=True
    )
    reference_id: Mapped[int] = mapped_column(
        db.ForeignKey("edit.references.id"), primary_key=True
    )
