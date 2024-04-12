from sqlalchemy.orm import Mapped, mapped_column

from submission.extensions import db
from submission.utils import ReferenceUtils


class Reference(db.Model):
    __tablename__ = "references"
    __table_args__ = {"schema": "edit"}

    reference_id: Mapped[int] = mapped_column(autoincrement=True, primary_key=True)
    doi: Mapped[str] = mapped_column(nullable=True)
    pubmed: Mapped[str] = mapped_column(nullable=True)
    title: Mapped[str] = mapped_column(nullable=True)
    authors: Mapped[str] = mapped_column(nullable=True)
    year: Mapped[int] = mapped_column(nullable=True)
    journal: Mapped[str] = mapped_column(nullable=True)

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

    @staticmethod
    def load(reference: str):
        """Loads a reference into the reference table

        Args:
            reference (str): reference formatted as doi:10... | pubmed:73...
        """
        metadata = ReferenceUtils.get_reference_metadata(reference)
        ref = Reference(
            doi=metadata.get("doi"),
            pubmed=metadata.get("pmid"),
            title=metadata.get("title"),
            authors=metadata.get("authors"),
            year=metadata.get("year"),
            journal=metadata.get("journal"),
        )
        db.session.add(ref)
        db.session.commit()
