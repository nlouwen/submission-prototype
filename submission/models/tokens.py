from sqlalchemy import insert, delete
from sqlalchemy.orm import Mapped, mapped_column
import datetime
import uuid

from submission.extensions import db


class Token(db.Model):
    __tablename__ = "tokens"

    token_id: Mapped[str] = mapped_column(
        primary_key=True, default=lambda: str(uuid.uuid4())
    )
    user_id: Mapped[int] = mapped_column(db.ForeignKey("users.id"), nullable=False)
    purpose: Mapped[str] = mapped_column(nullable=False)
    created_at: Mapped[datetime.datetime] = mapped_column(default=db.func.now())

    @staticmethod
    def generate_token(user_id: int, purpose: str) -> str:
        """Generate a token for a user with a specified purpose

        Args:
            user_id (int): db id of user
            purpose (str): desciption of token purpose
        """
        statement = (
            insert(Token)
            .values(user_id=user_id, purpose=purpose)
            .returning(Token.token_id)
        )
        ret = db.session.execute(statement).fetchone()

        if ret is None:
            raise RuntimeError("Database did not return token")

        db.session.commit()

        return ret[0]

    def cleanup_tokens(self) -> None:
        """Delete all tokens matching the user_id and purpose"""
        statement = delete(Token).where(
            Token.user_id == self.user_id, Token.purpose == self.purpose
        )
        db.session.execute(statement)
        db.session.commit()

    def is_created_within(self, hours: float) -> bool:
        """Check whether token was created within a number of hours

        Args:
            hours (float): number of hours

        Returns:
            bool: true if token was created within specified number of hours
        """
        current_time = datetime.datetime.now(datetime.UTC)
        boundary = current_time - datetime.timedelta(hours=hours)
        return self.created_at > boundary
