from wtforms import Form, SubmitField


class EditSelectForm(Form):
    minimal = SubmitField("Minimal Entry")
    structure = SubmitField("Structure")
    activity = SubmitField("Biological activity")
    biosynth = SubmitField("Class-specific and biosynthesis information")
    tailoring = SubmitField("Tailoring reaction information")
    annotation = SubmitField("Gene annotation")
