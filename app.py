from flask import (
    Flask,
    render_template,
    render_template_string,
    request,
    redirect,
    url_for,
    make_response,
    flash,
)
from forms.existing_bgc import SelectExisting
from forms.edit_select import EditSelectForm
from forms.common import is_valid_bgc_id
from forms.form_collection import FormCollection
from rdkit import Chem, rdBase
from rdkit.Chem import Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from werkzeug.datastructures import MultiDict
from pathlib import Path
import json
import csv
import time

app = Flask(__name__)
app.secret_key = "IYKYK"


### Landing page
@app.route("/", methods=["GET", "POST"])
def index():
    """Main page, edit existing entry or start new entry"""
    form = SelectExisting(request.form)
    if request.method == "POST":
        # create new entry
        if form.submit.data:
            bgc_id = create_new_entry()
            return redirect(url_for("edit_bgc", bgc_id=bgc_id))

        # edit valid existing entry
        if request.form.get("edit") and form.validate():
            bgc_id = form.accession.data
            Path(f"{bgc_id}_data.json").touch()

            return redirect(url_for("edit_bgc", bgc_id=bgc_id))

    return render_template("index.html", form=form)


## TODO: Submit new
@app.route("/WIP")
def work_in_progress():
    return render_template("WIP.html")


## Edit existing
@app.route("/edit/<bgc_id>", methods=["GET", "POST"])
def edit_bgc(bgc_id: str):

    if not is_valid_bgc_id(bgc_id):
        return "Invalid existing entry!", 404

    form = EditSelectForm(request.form)
    if request.method == "POST":
        for option, chosen in form.data.items():
            if chosen:
                return redirect(url_for(f"edit_{option}", bgc_id=bgc_id))
    return render_template("edit.html", form=form, bgc_id=bgc_id)


# Edit minimal information
@app.route("/edit/<bgc_id>/minimal", methods=["GET", "POST"])
def edit_minimal(bgc_id: str):
    if not is_valid_bgc_id(bgc_id):
        return "Invalid existing entry!", 404

    # try to fill data from existing entry
    if not request.form:
        form = FormCollection.minimal(MultiDict(read_data(bgc_id).get("Minimal")))
    else:
        form = FormCollection.minimal(request.form)

    if request.method == "POST" and form.validate():
        # TODO: save to database
        save_data(bgc_id, "Minimal", request.form)
        flash("Submitted minimal entry!")
        return redirect(url_for("edit_bgc", bgc_id=bgc_id))
    return render_template("min_entry.html", form=form, bgc_id=bgc_id)


# Edit structure information
@app.route("/edit/<bgc_id>/structure", methods=["GET", "POST"])
def edit_structure(bgc_id: str):
    """Main page for editing entry structure information"""
    if not is_valid_bgc_id(bgc_id):
        return "Invalid existing entry!", 404

    if not request.form:
        form = FormCollection.structure(MultiDict(read_data(bgc_id).get("Structure")))
    else:
        form = FormCollection.structure(request.form)
    if request.method == "POST":
        if form.add.data:
            form.structures.append_entry()
            return render_template("structure.html", form=form, bgc_id=bgc_id)
        elif form.submit.data and form.validate():
            # return form.data  # TODO: save to db
            save_data(bgc_id, "Structure", request.form)
            flash("Submitted structure information!")
            return redirect(url_for("edit_bgc", bgc_id=bgc_id))

    # on GET query db for any products already present
    else:
        # if one (or more) is present prefill
        try:
            products = next(
                csv.reader(
                    [MultiDict(read_data(bgc_id).get("Minimal")).get("products")],
                    skipinitialspace=True,
                )
            )
        except:
            products = [""]

        for product in products:
            if product not in [struct.data.get("name") for struct in form.structures]:
                form.structures.append_entry(data={"name": product})
    return render_template("structure.html", form=form, bgc_id=bgc_id)


@app.route("/render-smiles", methods=["POST"])
def render_smiles():
    """Render entered smiles string using rdkit"""
    # record which product smiles was entered
    origin = request.headers["Hx-Trigger-Name"]
    smiles = request.form.get(origin).strip()

    if not smiles:
        return ""

    # RDkit implementation
    try:
        mol = Chem.MolFromSmiles(smiles)
        rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 200)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace("svg:", "")
    except Exception:
        return "Invalid SMILES"

    response = make_response(svg)
    response.content_type = "image/svg+xml"
    return response


## Bio Act
@app.route("/edit/<bgc_id>/bioact", methods=["GET", "POST"])
def edit_activity(bgc_id: str):
    if not is_valid_bgc_id(bgc_id):
        return "Invalid existing entry!", 404

    if not request.form:
        form = FormCollection.bioact(MultiDict(read_data(bgc_id).get("Bio_activity")))
    else:
        form = FormCollection.bioact(request.form)
    if request.method == "POST":
        # for activity in form.activities:
        #     if activity.add.data:
        #         activity.assay.append_entry()
        #         return render_template(
        #             "biological_activity.html", bgc_id=bgc_id, form=form
        #         )
        if form.submit.data and form.validate():
            # return form.data # TODO: save to db
            save_data(bgc_id, "Bio_activity", request.form)
            flash("Submitted activity information!")
            return redirect(url_for("edit_bgc", bgc_id=bgc_id))
    else:
        # prefill compounds
        try:
            products = next(
                csv.reader(
                    [MultiDict(read_data(bgc_id).get("Minimal")).get("products")],
                    skipinitialspace=True,
                )
            )
        except:
            products = [""]

        for product in products:
            if product not in [act.data.get("compound") for act in form.activities]:
                form.activities.append_entry(data={"compound": product})
    return render_template("biological_activity.html", bgc_id=bgc_id, form=form)


## Biosynthesis
@app.route("/edit/<bgc_id>/biosynth", methods=["GET", "POST"])
def edit_biosynth(bgc_id: str):
    if not is_valid_bgc_id(bgc_id):
        return "Invalid existing entry!", 404

    # if request.method == "POST":
    #     b_class = request.form.get("b_class")
    #     form = getattr(FormCollection, b_class)(request.form)

    #     return form.data
    # return request.form

    # form = BiosyntheticClassesForm(request.form)
    return render_template("biosynthesis.html", bgc_id=bgc_id)


@app.route("/edit/<bgc_id>/biosynth/<b_class>", methods=["GET", "POST"])
def edit_biosynth_class(bgc_id: str, b_class: str):
    if not is_valid_bgc_id(bgc_id):
        return "Invalid existing entry!", 404

    if not request.form:
        form = getattr(FormCollection, b_class)(
            MultiDict(read_data(bgc_id).get(f"BioSynth_{b_class}"))
        )
    else:
        form = getattr(FormCollection, b_class)(request.form)

    if request.method == "POST":

        # Add extra NRPS fields
        # if form.data.get("add_release_type"):
        #     form.release_types.append_entry()
        #     return render_template(
        #         "biosynth_class_specific.html",
        #         form=form,
        #         b_class=b_class,
        #         bgc_id=bgc_id,
        #     )
        # if form.data.get("add_thioesterase"):
        #     form.thioesterases.append_entry()
        #     return render_template(
        #         "biosynth_class_specific.html",
        #         form=form,
        #         b_class=b_class,
        #         bgc_id=bgc_id,
        #     )

        # Add extra Ribosomal fields
        # if form.data.get("precursors"):
        #     if form.data.get("add_precursors"):
        #         form.precursors.append_entry()
        #         return render_template(
        #             "biosynth_class_specific.html",
        #             form=form,
        #             b_class=b_class,
        #             bgc_id=bgc_id,
        #         )
        #     for precursor in form.precursors:
        #         if precursor.data.get("add_crosslinks"):
        #             precursor.crosslinks.append_entry()
        #             return render_template(
        #                 "biosynth_class_specific.html",
        #                 form=form,
        #                 b_class=b_class,
        #                 bgc_id=bgc_id,
        #             )

        # # Add extra Saccharide fields
        # if form.data.get("add_glycosyltransferase"):
        #     form.glycosyltransferases.append_entry()
        #     return render_template(
        #         "biosynth_class_specific.html",
        #         form=form,
        #         b_class=b_class,
        #         bgc_id=bgc_id,
        #     )
        # if form.data.get("add_subcluster"):
        #     form.subclusters.append_entry()
        #     return render_template(
        #         "biosynth_class_specific.html",
        #         form=form,
        #         b_class=b_class,
        #         bgc_id=bgc_id,
        #     )

        # if no adds were triggered, validate and process data
        if form.validate():
            save_data(bgc_id, f"BioSynth_{b_class}", request.form)
            flash(f"Submitted {b_class} biosynthesis information!")
            return redirect(url_for("edit_bgc", bgc_id=bgc_id))
            # return form.data

    return render_template(
        "biosynth_class_specific.html", form=form, b_class=b_class, bgc_id=bgc_id
    )


# @app.route("/get_class", methods=["POST"])
# def get_class():
#     if request.headers.get("Hx-Request"):
#         b_class = request.form.get("b_class")
#         form = getattr(BioClassesCollection, b_class)()
#         return render_template("forms_basic.html", form=form)


@app.route("/class_buttons/<bgc_id>", methods=["POST"])
def class_buttons(bgc_id: str):

    # grab classes for current bgc_id
    classes = MultiDict(read_data(bgc_id).get("Minimal")).getlist("b_class")

    if not classes:
        classes = ["NRPS", "PKS", "Ribosomal", "Saccharide", "Terpene", "Other"]
    class_btns = ""
    for cls in classes:
        # class_btns += f"<a href='/edit/{bgc_id}/biosynth/{cls}'><button class='btn btn-light'>{cls}</button><a>"
        class_btns += f"<a class='btn btn-light' style='margin: 5px' role='button' href='/edit/{bgc_id}/biosynth/{cls}'>{cls}</a>"
    return class_btns


## tailoring
@app.route("/edit/<bgc_id>/tailoring", methods=["GET", "POST"])
def edit_tailoring(bgc_id: str):
    if not is_valid_bgc_id(bgc_id):
        return "Invalid existing entry!", 404

    if not request.form:
        form = FormCollection.tailoring(MultiDict(read_data(bgc_id).get("Tailoring")))
    else:
        form = FormCollection.tailoring(request.form)

    if request.method == "POST" and form.validate():
        save_data(bgc_id, "Tailoring", request.form)
        flash("Submitted tailoring information!")
        return redirect(url_for("edit_bgc", bgc_id=bgc_id))

    return render_template("tailoring.html", bgc_id=bgc_id, form=form)


@app.route("/edit/<bgc_id>/annotation", methods=["GET", "POST"])
def edit_annotation(bgc_id: str):
    if not is_valid_bgc_id(bgc_id):
        return "Invalid existing entry!", 404

    if not request.form:
        form = FormCollection.annotation(MultiDict(read_data(bgc_id).get("Annotation")))
    else:
        form = FormCollection.annotation(request.form)

    if request.method == "POST" and form.validate():
        save_data(bgc_id, "Annotation", request.form)
        flash("Submitted annotation information!")
        return redirect(url_for("edit_bgc", bgc_id=bgc_id))

    return render_template("annotation.html", bgc_id=bgc_id, form=form)


## utils
# TODO: generalize to obtain present genes etc.
@app.route("/query_ncbi", methods=["POST"])
def query_ncbi():
    """Very rough outline of ncbi query to obtain taxonomy info"""
    # mock flow
    time.sleep(5)

    try:
        # accession = request.form.get("genome")
        # query ncbi
        # ncbi_record = Entrez.read(
        #                 Entrez.efetch(db="nuccore", id=accession, retmode="xml")
        #             )
        # organism = ncbi_record[0][
        #                 "GBSeq_source"
        #             ]
        organism = "Streptomyces coelicolor A3(2)"

        # naively parse
        genus, species, *strain = organism.split()

        # look up organism in ncbi names.dmp for tax id
        tax_id = 100226
        message = "Autofilled based on accession, please doublecheck that all data was filled correctly."
    except:
        # upon error encountered
        tax_id, genus, species, strain = ("", "", "", [])
        message = "Unable to find taxonomy information based on accession, please enter manually."

    formdata = MultiDict(
        [
            ("ncbi_tax_id", tax_id),
            ("genus", genus),
            ("species", species),
            ("strain", " ".join(strain)),
        ]
    )
    form = FormCollection.minimal.TaxonomyForm(formdata)
    return render_template_string(
        """{% import 'macros.html' as m %}
        <span class="form-text text-muted fst-italic">{{message}}</span>
        {{m.simple_divsubform(field)}}""",
        field=form,
        message=message,
    )


# minimal
@app.route("/add_locus", methods=["POST"])
def add_locus():
    form = FormCollection.minimal(request.form)
    form.loci.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=form.loci[-1],
    )


@app.route("/add_evidence", methods=["POST"])
def add_evidence():
    form = FormCollection.minimal(request.form)
    _, locus_idx, _ = request.headers.get("Hx-Trigger").split("-")
    locus = form.loci[int(locus_idx)]
    locus.evidence.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=locus.evidence[-1],
    )


# structures
@app.route("/add_compound", methods=["POST"])
def add_compound():
    form = FormCollection.structure(request.form)
    form.structures.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, struct=True, deletebtn=True)}}""",
        field=form.structures[-1],
    )


# Bio activities
@app.route("/add_assay", methods=["POST"])
def add_assay():
    form = FormCollection.bioact(request.form)
    _, origin_idx, _ = request.headers.get("Hx-Trigger").split("-")
    origin_activity = form.activities[int(origin_idx)]
    origin_activity.assays.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=origin_activity.assays[-1],
    )


@app.route("/add_bioact_compound", methods=["POST"])
def add_bioact_compound():
    form = FormCollection.bioact(request.form)
    form.activities.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=form.activities[-1],
    )


# NRPS
@app.route("/add_release", methods=["POST"])
def add_release():
    form = FormCollection.NRPS(request.form)
    form.release_types.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=form.release_types[-1],
    )


@app.route("/add_thioesterase", methods=["POST"])
def add_thioesterase():
    form = FormCollection.NRPS(request.form)
    form.thioesterases.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=form.thioesterases[-1],
    )


# ribosomal
@app.route("/add_precursor", methods=["POST"])
def add_precursor():
    form = FormCollection.Ribosomal(request.form)
    form.precursors.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=form.precursors[-1],
    )


@app.route("/add_crosslink", methods=["POST"])
def add_crosslink():
    form = FormCollection.Ribosomal(request.form)
    _, origin_idx, _ = request.headers.get("Hx-Trigger").split("-")
    origin_precursor = form.precursors[int(origin_idx)]
    origin_precursor.crosslinks.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=origin_precursor.crosslinks[-1],
    )


# saccharide
@app.route("/add_glycosyltransferase", methods=["POST"])
def add_glycosyltransferase():
    form = FormCollection.Saccharide(request.form)
    form.glycosyltransferases.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=form.glycosyltransferases[-1],
    )


@app.route("/add_subcluster", methods=["POST"])
def add_subcluster():
    form = FormCollection.Saccharide(request.form)
    form.subclusters.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=form.subclusters[-1],
    )


# tailoring
@app.route("/add_tailoring_enzyme", methods=["POST"])
def add_tailoring_enzyme():
    form = FormCollection.tailoring(request.form)
    form.enzymes.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=form.enzymes[-1],
    )


@app.route("/add_tailoring_reaction", methods=["POST"])
def add_tailoring_reaction():
    form = FormCollection.tailoring(request.form)
    _, origin_idx, _ = request.headers.get("Hx-Trigger").split("-")
    origin_enzyme = form.enzymes[int(origin_idx)]
    origin_enzyme.reactions.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=origin_enzyme.reactions[-1],
    )


@app.route("/add_aux_enzyme", methods=["POST"])
def add_aux_enzyme():
    form = FormCollection.tailoring(request.form)
    _, origin_idx, _ = request.headers.get("Hx-Trigger").split("-", 2)
    origin_enzyme = form.enzymes[int(origin_idx)]
    origin_enzyme.enzyme.auxiliary_enzymes.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=origin_enzyme.enzyme.auxiliary_enzymes[-1],
    )


@app.route("/add_ontology", methods=["POST"])
def add_ontology():
    form = FormCollection.tailoring(request.form)
    _, enzyme_idx, _, reaction_idx, _ = request.headers.get("Hx-Trigger").split("-")
    enzyme = form.enzymes[int(enzyme_idx)]
    reaction = enzyme.reactions[int(reaction_idx)]
    reaction.tailoring.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=reaction.tailoring[-1],
    )


@app.route("/add_smarts", methods=["POST"])
def add_smarts():
    form = FormCollection.tailoring(request.form)
    _, enzyme_idx, _, reaction_idx, _ = request.headers.get("Hx-Trigger").split("-")
    enzyme = form.enzymes[int(enzyme_idx)]
    reaction = enzyme.reactions[int(reaction_idx)]
    reaction.reaction_smarts.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=reaction.reaction_smarts[-1],
    )


@app.route("/add_val_reaction", methods=["POST"])
def add_val_reaction():
    form = FormCollection.tailoring(request.form)
    _, enzyme_idx, _, reaction_idx, _ = request.headers.get("Hx-Trigger").split("-")
    enzyme = form.enzymes[int(enzyme_idx)]
    reaction = enzyme.reactions[int(reaction_idx)]
    reaction.validated_reactions.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=reaction.validated_reactions[-1],
    )


@app.route("/add_hydrogen", methods=["POST"])
def add_hydrogen():
    form = FormCollection.tailoring(request.form)
    _, enzyme_idx, _, reaction_idx, _, smarts_idx, _ = request.headers.get(
        "Hx-Trigger"
    ).split("-")
    enzyme = form.enzymes[int(enzyme_idx)]
    reaction = enzyme.reactions[int(reaction_idx)]
    smarts = reaction.reaction_smarts[int(smarts_idx)]
    smarts.explicitHydrogen.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=smarts.explicitHydrogen[-1],
    )


@app.route("/add_tail_reaction_smarts_evidence", methods=["POST"])
def add_tail_reaction_smarts_evidence():
    form = FormCollection.tailoring(request.form)
    _, enzyme_idx, _, reaction_idx, _, smarts_idx, _ = request.headers.get(
        "Hx-Trigger"
    ).split("-")
    enzyme = form.enzymes[int(enzyme_idx)]
    reaction = enzyme.reactions[int(reaction_idx)]
    smarts = reaction.reaction_smarts[int(smarts_idx)]
    smarts.evidence_sm.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=smarts.evidence_sm[-1],
    )


@app.route("/add_val_reaction_evidence", methods=["POST"])
def add_val_reaction_evidence():
    form = FormCollection.tailoring(request.form)
    _, enzyme_idx, _, reaction_idx, _, val_idx, _ = request.headers.get(
        "Hx-Trigger"
    ).split("-")
    enzyme = form.enzymes[int(enzyme_idx)]
    reaction = enzyme.reactions[int(reaction_idx)]
    val = reaction.reaction_smarts[int(val_idx)]
    val.evidence_val.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=val.evidence_val[-1],
    )


@app.route("/add_general", methods=["POST"])
def add_general():
    directions = request.headers["Hx-Trigger"].split("-")
    formname = Path(request.referrer).name

    curr = getattr(FormCollection, formname)(request.form)
    for i in range(0, len(directions), 2):
        if i + 2 > len(directions):
            final = getattr(curr, directions[i])
            final.append_entry()
            return render_template_string(
                """{% import 'macros.html' as m %}
                {{m.simple_divsubform(field, deletebtn=True)}}""",
                field=final[-1],
            )

        else:
            subform, subform_idx = directions[i : i + 2]
            curr = getattr(curr, subform)[int(subform_idx)]


@app.route("/delete", methods=["DELETE"])
def delete():
    return ""


@app.route("/submit", methods=["GET", "POST"])
def submit():
    return render_template("base.html")


## Temp save data to file
def save_data(bgc_id: str, section_key: str, req_data: MultiDict):
    """Append data to file containing all answers for one BGC"""
    data = {section_key: [(k, v) for k in req_data for v in req_data.getlist(k)]}
    existing_data = read_data(bgc_id)

    existing_data.update(data)
    with open(f"{bgc_id}_data.json", "w") as outf:
        json.dump(existing_data, outf, sort_keys=True, indent=4)


def read_data(bgc_id: str):
    with open(f"{bgc_id}_data.json", "r") as inf:
        if content := inf.read():
            existing_data = json.loads(content)
        else:
            existing_data = {}
    return existing_data


def create_new_entry():
    max_entry_id = 0
    for filepath in Path(".").glob("new*_data.json"):
        if (nr := int(filepath.stem[3:6])) > max_entry_id:
            max_entry_id = nr
    bgc_id = f"new{max_entry_id+1:0>3}"
    Path(f"{bgc_id}_data.json").touch()

    return bgc_id
