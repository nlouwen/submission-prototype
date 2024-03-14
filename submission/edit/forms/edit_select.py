from wtforms import Form, SubmitField


class EditSelectForm(Form):
    minimal = SubmitField("Minimal Entry")
    structure = SubmitField("Structure")
    activity = SubmitField("Biological activity")
    biosynth = SubmitField("Class-specific information")
    tailoring = SubmitField("Tailoring reaction information")
    annotation = SubmitField("Gene annotation")
