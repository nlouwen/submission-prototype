"""empty message

Revision ID: affa9a9fb156
Revises: 23ed6242ea5e
Create Date: 2024-03-15 13:41:28.484213

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'affa9a9fb156'
down_revision = '23ed6242ea5e'
branch_labels = None
depends_on = None


def upgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.create_table('user_info',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('alias', sa.String(), nullable=False),
    sa.Column('name', sa.String(), nullable=False),
    sa.Column('call_name', sa.String(), nullable=False),
    sa.Column('organisation', sa.String(), nullable=False),
    sa.Column('public', sa.Boolean(), nullable=False),
    sa.ForeignKeyConstraint(['id'], ['users.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('alias')
    )
    with op.batch_alter_table('users', schema=None) as batch_op:
        batch_op.add_column(sa.Column('active', sa.Boolean(), nullable=False))
        batch_op.drop_column('name')

    # ### end Alembic commands ###


def downgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    with op.batch_alter_table('users', schema=None) as batch_op:
        batch_op.add_column(sa.Column('name', sa.VARCHAR(), nullable=False))
        batch_op.drop_column('active')

    op.drop_table('user_info')
    # ### end Alembic commands ###
