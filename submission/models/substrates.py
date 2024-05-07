from sqlalchemy import select, inspect, or_, ScalarResult
from sqlalchemy.orm import Mapped, mapped_column

from submission.extensions import db


class Substrate(db.Model):
    """Static collection of substrates"""

    __tablename__ = "substrates"
    __table_args__ = {"schema": "edit"}

    id: Mapped[int] = mapped_column(autoincrement=True, primary_key=True)
    name: Mapped[str]
    fullname: Mapped[str] = mapped_column(nullable=True)
    structure: Mapped[str]
    proteinogenic: Mapped[bool]

    @property
    def identifier(self):
        if not self.fullname:
            return self.name
        return self.fullname

    def summarize(self, html=True) -> str:
        if self.fullname:
            full = f": {self.fullname}"
        else:
            full = ""
        if html:
            return f"<b>{self.name}</b>{full}"
        else:
            return f"{self.name}{full}"

    def __repr__(self) -> str:
        return self.summarize(html=False)

    @staticmethod
    def isearch(value) -> ScalarResult["Substrate"]:
        """Search case insensitive with matches anywhere in the name or fullname

        Args:
            value (str): search term

        Returns:
            ScalarResult[Substrate]: collection of substrates matching search term, can be empty
        """
        return db.session.scalars(
            select(Substrate).where(
                or_(
                    Substrate.name.icontains(value), Substrate.fullname.icontains(value)
                )
            )
        )
