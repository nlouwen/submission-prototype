from typing import Literal, Optional, Union

import requests
from sqlalchemy import select
from sqlalchemy.orm import Mapped, mapped_column

from submission.extensions import db
from submission.utils.custom_errors import ReferenceNotFound


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
    url: Mapped[str] = mapped_column(nullable=True)
    patent: Mapped[str] = mapped_column(nullable=True)
    entries = db.relationship(
        "Entry",
        secondary="edit.entry_references",
        back_populates="references",
        lazy="selectin",
    )

    @property
    def identifier(self):
        """Access a reference identifier, prefers doi>pubmed>url>patent respectively

        If somehow no identifier is present, returns None
        """
        if self.doi:
            return f"doi:{self.doi}"
        elif self.pubmed:
            return f"pubmed:{self.pubmed}"
        elif self.url:
            return f"url:{self.url}"
        elif self.patent:
            return f"patent:{self.patent}"
        else:
            return None

    def __repr__(self):
        return self.summarize(html=False)

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

    def summarize(self, html=True):
        """Generate a one-line summary of the reference"""
        if self.doi == "pending":
            return "Pending Publication Placeholder"
        if self.url or self.patent:
            return self.identifier
        title = self.title or ""
        if not title.endswith("."):
            title += "."
        journal = self.journal or ""
        year = self.year or ""
        identifier = self.identifier or ""
        if html:
            return f"{title} {self.short_authors()} <i>{journal}</i>, <b>{year}</b>. {identifier}"
        return f"{title} {self.short_authors()} {journal}, {year}. {identifier}"

    @classmethod
    def load(cls, reference: str):
        """Loads a reference into the reference table

        Args:
            reference (str): reference formatted as doi:10... | pubmed:73...
        """
        if reference.startswith("doi:") or reference.startswith("pubmed:"):
            # distinguish Biorxiv and Chemrxiv dois
            if reference.startswith("doi:10.1101/"):
                metadata = Reference.get_biorxiv_metadata(reference)
            elif reference.startswith("doi:10.26434/"):
                metadata = Reference.get_chemrxiv_metadata(reference)
            else:
                metadata = Reference.get_reference_metadata(reference)

            if not metadata or metadata.get("detail") == "Article not found":
                raise ReferenceNotFound(reference)

            ref = cls(
                doi=metadata.get("doi"),
                pubmed=metadata.get("pmid"),
                title=metadata.get("title"),
                authors=metadata.get("authors"),
                year=metadata.get("year"),
                journal=metadata.get("journal"),
            )
        elif reference.startswith("url"):
            ref = cls(url=reference.split(":", 1)[1])
        elif reference.startswith("patent"):
            ref = cls(patent=reference.split(":", 1)[1])
        else:
            raise ReferenceNotFound(reference)
        db.session.add(ref)
        db.session.commit()
        return ref

    @staticmethod
    def load_missing(references: list[str]) -> list["Reference"]:
        """Load missing references into database

        Args:
            references (list[str]): references to load, fmt: 'doi:10..' | 'pubmed:77..'

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
        id_type, ident = reference.split(":", 1)
        if id_type == "doi":
            return db.session.scalar(select(Reference).where(Reference.doi == ident))
        elif id_type == "pubmed":
            return db.session.scalar(select(Reference).where(Reference.pubmed == ident))
        elif id_type == "url":
            return db.session.scalar(select(Reference).where(Reference.url == ident))
        elif id_type == "patent":
            return db.session.scalar(select(Reference).where(Reference.patent == ident))
        else:
            raise RuntimeError(f"Unexpected reference type '{id_type}'")

    @staticmethod
    def get_reference_metadata(reference: str) -> dict[str, Optional[str]]:
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
        if r.status_code != 200:
            return {}
        return r.json()

    @staticmethod
    def get_biorxiv_metadata(reference: str) -> dict[str, Optional[str]]:
        """Collect metadata for Biorxiv references using Biorxiv API

        Args:
            reference (str): reference formatted as doi:10.1101/...

        Returns:
            dict[str, str]: mapping of reference metadata
        """

        def post_biorxiv_server(
            server: Literal["biorxiv", "medrxiv"], identifier: str
        ) -> dict[str, Optional[str]]:
            """POST request to biorxiv api

            Args:
                server (str): server to query, biorxiv or medrxiv
                identifier (str): doi identifier, starts with 10.1101/

            Returns:
                dict[str, str]: mapping of reference metadata if reference was found
            """
            r = requests.post(f"https://api.biorxiv.org/details/{server}/{identifier}")
            if r.status_code != 200:
                return {}

            raw: dict = r.json()
            if raw["messages"][0]["status"] != "ok":
                return {}
            else:
                # take the latest version of reference
                mapping: dict = raw["collection"][-1]
                mapping["year"] = (
                    d.split("-", 1)[0] if (d := mapping.get("date")) else None
                )
                mapping["journal"] = mapping.get("server")
                return mapping

        identifier = reference.split(":", 1)[1]
        mapping = post_biorxiv_server("biorxiv", identifier)
        if not mapping:
            # try medrxiv to be sure
            mapping = post_biorxiv_server("medrxiv", identifier)
        return mapping

    @staticmethod
    def get_chemrxiv_metadata(reference: str) -> dict[str, Optional[str]]:
        """Collect metadata for Chemrxiv references using the Chemrxiv API

        Args:
            reference (str): reference formatted as doi:10.26434/...

        Returns:
            dict[str, str]: mapping of reference metadata
        """
        identifier = reference.split(":", 1)[1]
        r = requests.get(
            f"https://chemrxiv.org/engage/chemrxiv/public-api/v1/items/doi/{identifier}"
        )
        # not found returns 404 status code
        if r.status_code != 200:
            return {}

        raw: dict = r.json()
        if pub_d := raw.get("publishedDate"):
            year = pub_d.split("-")[0]
        else:
            year = None

        mapping = {
            "doi": raw.get("doi"),
            "title": raw.get("title"),
            "authors": "; ".join(
                f'{author.get("lastName")}, {author.get("firstName")}'
                for author in raw["authors"]
            ),
            "year": year,
            "journal": "chemrxiv",
        }
        return mapping


class EntryReference(db.Model):
    __tablename__ = "entry_references"
    __table_args__ = {"schema": "edit"}

    entry_id: Mapped[int] = mapped_column(
        db.ForeignKey("edit.entries.id"), primary_key=True
    )
    reference_id: Mapped[int] = mapped_column(
        db.ForeignKey("edit.references.id"), primary_key=True
    )
