from flask_login import UserMixin
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import Mapped, mapped_column
from werkzeug.security import generate_password_hash, check_password_hash

from extensions import db


class User(UserMixin, db.Model):
    __tablename__ = "users"
    #__table_args__ = {"schema": "auth"}

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(nullable=False)
    email: Mapped[str] = mapped_column(unique=True, nullable=False)
    _password: Mapped[str] = mapped_column(nullable=False)

    roles: Mapped[list["Role"]] = db.relationship("Role", secondary="user_roles", back_populates="users", lazy="selectin")

    @hybrid_property
    def password(self):
        return self._password

    @password.setter
    def password(self, password):
        self._password = generate_password_hash(password)

    def check_password(self, password) -> bool:
        return check_password_hash(self.password, password)

    def has_role(self, role) -> bool:
        return bool(
            Role.query
            .join(Role.users)
            .filter(User.id == self.id)
            .filter(Role.slug == role)
            .count() == 1
        )


class Role(db.Model):
    __tablename__ = "roles"
    #__table_args__ = {"schema": "auth"}

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(nullable=False)
    slug: Mapped[str] = mapped_column(nullable=False, unique=True)

    users: Mapped[list["User"]] = db.relationship("User", secondary="user_roles", back_populates="roles", lazy="selectin")

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return f"Role('{self.name}', '{self.slug}')"


class UserRole(db.Model):
    __tablename__ = "user_roles"
    #__table_args__ = {"schema": "auth"}

    user_id: Mapped[int] = mapped_column(db.ForeignKey("users.id"), primary_key=True)
    role_id: Mapped[int] = mapped_column(db.ForeignKey("roles.id"), primary_key=True)
