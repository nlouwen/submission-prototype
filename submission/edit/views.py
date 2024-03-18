import re
import csv
from pathlib import Path

from flask import (
    abort,
    render_template,
    render_template_string,
    request,
    redirect,
    url_for,
    flash,
)
from flask_login import current_user, login_required
from werkzeug.datastructures import MultiDict
from werkzeug.wrappers import response
from wtforms.widgets import html_params
from markupsafe import Markup

from submission.edit import bp_edit
from submission.edit.forms.form_collection import FormCollection
from submission.edit.forms.edit_select import EditSelectForm
from submission.utils import Storage, draw_smiles_svg, ReferenceUtils
from submission.utils.custom_validators import is_valid_bgc_id


@bp_edit.route("/<bgc_id>", methods=["GET", "POST"])
@login_required
def edit_bgc(bgc_id: str) -> str | response.Response:
    """Overview page with navigation to forms for entry sections

    Args:
        bgc_id (str): BGC identifier

    Returns:
        str | Response: rendered overview template or redirect section form
    """

    if not is_valid_bgc_id(bgc_id):
        return abort(403, "Invalid existing entry!")

    form = EditSelectForm(request.form)
    if request.method == "POST":
        for option, chosen in form.data.items():
            if chosen:
                return redirect(url_for(f"edit.edit_{option}", bgc_id=bgc_id))
    return render_template("edit/edit.html", form=form, bgc_id=bgc_id)


@bp_edit.route("/<bgc_id>/minimal", methods=["GET", "POST"])
@login_required
def edit_minimal(bgc_id: str) -> str | response.Response:
    """Form to enter minimal entry information

    Args:
        bgc_id (str): BGC identifier

    Returns:
        str | Response: rendered form template or redirect to edit_bgc overview
    """
    if not is_valid_bgc_id(bgc_id):
        return abort(403, "Invalid existing entry!")

    # try to fill data from existing entry
    if not request.form:
        form = FormCollection.minimal(
            MultiDict(Storage.read_data(bgc_id).get("Minimal"))
        )
    else:
        form = FormCollection.minimal(request.form)

    if request.method == "POST" and form.validate():
        # TODO: save to database
        Storage.save_data(bgc_id, "Minimal", request.form)
        flash("Submitted minimal entry!")
        return redirect(url_for("edit.edit_bgc", bgc_id=bgc_id))
    return render_template(
        "edit/min_entry.html",
        form=form,
        bgc_id=bgc_id,
        is_reviewer=current_user.has_role("reviewer"),
    )


@bp_edit.route("/<bgc_id>/structure", methods=["GET", "POST"])
@login_required
def edit_structure(bgc_id: str) -> str | response.Response:
    """Form to enter structure information

    Args:
        bgc_id (str): BGC identifier

    Returns:
        str | Response: rendered form template or redirect to edit_bgc overview
    """
    if not is_valid_bgc_id(bgc_id):
        return abort(403, "Invalid existing entry!")

    if not request.form:
        form = FormCollection.structure(
            MultiDict(Storage.read_data(bgc_id).get("Structure"))
        )
    else:
        form = FormCollection.structure(request.form)

    if request.method == "POST" and form.validate():
        # TODO: save to db
        Storage.save_data(bgc_id, "Structure", request.form)
        flash("Submitted structure information!")
        return redirect(url_for("edit.edit_bgc", bgc_id=bgc_id))

    # on GET query db for any products already present
    else:
        # if one (or more) is present prefill
        try:
            products = next(
                csv.reader(
                    [
                        MultiDict(Storage.read_data(bgc_id).get("Minimal")).get(
                            "products"
                        )
                    ],
                    skipinitialspace=True,
                )
            )
        except:
            products = []

        for product in products:
            if product not in [struct.data.get("name") for struct in form.structures]:
                form.structures.append_entry(data={"name": product})
    return render_template(
        "edit/structure.html",
        form=form,
        bgc_id=bgc_id,
        is_reviewer=current_user.has_role("reviewer"),
    )


@bp_edit.route("/render_smiles", methods=["POST"])
def render_smiles() -> str | response.Response:
    origin = request.headers["Hx-Trigger-Name"]
    smiles_string = request.form.get(origin)

    if smiles_string is None or not (smiles := smiles_string.strip()):
        return ""

    return draw_smiles_svg(smiles)


@bp_edit.route("/<bgc_id>/bioact", methods=["GET", "POST"])
@login_required
def edit_activity(bgc_id: str) -> str | response.Response:
    """Form to enter biological activity information

    Args:
        bgc_id (str): BGC identifier

    Returns:
        str | Response: rendered form templare or redirect to edit_bgc overview
    """
    if not is_valid_bgc_id(bgc_id):
        return abort(403, "Invalid existing entry!")

    if not request.form:
        form = FormCollection.bioact(
            MultiDict(Storage.read_data(bgc_id).get("Bio_activity"))
        )
    else:
        form = FormCollection.bioact(request.form)
    if request.method == "POST" and form.validate():
        # TODO: save to db
        Storage.save_data(bgc_id, "Bio_activity", request.form)
        flash("Submitted activity information!")
        return redirect(url_for("edit.edit_bgc", bgc_id=bgc_id))

    else:
        # prefill compounds
        try:
            products = next(
                csv.reader(
                    [
                        MultiDict(Storage.read_data(bgc_id).get("Minimal")).get(
                            "products"
                        )
                    ],
                    skipinitialspace=True,
                )
            )
        except:
            products = []

        for product in products:
            if product not in [act.data.get("compound") for act in form.activities]:
                form.activities.append_entry(data={"compound": product})
    return render_template(
        "edit/biological_activity.html",
        bgc_id=bgc_id,
        form=form,
        is_reviewer=current_user.has_role("reviewer"),
    )


@bp_edit.route("/<bgc_id>/biosynth", methods=["GET", "POST"])
@login_required
def edit_biosynth(bgc_id: str) -> str:
    """Selection overview page for class-specific biosynthesis forms

    Args:
        bgc_id (str): BGC identifier

    Returns:
        str: rendered template
    """
    if not is_valid_bgc_id(bgc_id):
        return abort(403, "Invalid existing entry!")

    return render_template("edit/biosynthesis.html", bgc_id=bgc_id)


@bp_edit.route("/class_buttons/<bgc_id>", methods=["POST"])
def class_buttons(bgc_id: str) -> str:
    """Obtain buttons linking to relevant biosynthetic classes for this BGC

    Args:
        bgc_id (str): BGC identifier

    Returns:
        str: html buttons linking to class-specific form templates
    """

    # grab classes for current bgc_id
    classes = MultiDict(Storage.read_data(bgc_id).get("Minimal")).getlist("b_class")

    if not classes:
        classes = ["NRPS", "PKS", "Ribosomal", "Saccharide", "Terpene", "Other"]
    class_btns = ""
    for cls in classes:
        class_btns += f"<a class='btn btn-light' style='margin: 5px' role='button' href='/edit/{bgc_id}/biosynth/{cls}'>{cls}</a>"
    return class_btns


@bp_edit.route("/<bgc_id>/biosynth/<b_class>", methods=["GET", "POST"])
@login_required
def edit_biosynth_class(bgc_id: str, b_class: str) -> str | response.Response:
    """Form to enter class-specific biosynthesis information

    Args:
        bgc_id (str): BGC identifier
        b_class (str): Biosynthetic class

    Returns:
        str | Response: rendered template or redirect to edit_bgc overview
    """
    if not is_valid_bgc_id(bgc_id):
        return abort(403, "Invalid existing entry!")

    if not request.form:
        form = getattr(FormCollection, b_class)(
            MultiDict(Storage.read_data(bgc_id).get(f"BioSynth_{b_class}"))
        )
    else:
        form = getattr(FormCollection, b_class)(request.form)

    if request.method == "POST" and form.validate():
        # TODO: save to db
        Storage.save_data(bgc_id, f"BioSynth_{b_class}", request.form)
        flash(f"Submitted {b_class} biosynthesis information!")
        return redirect(url_for("edit.edit_bgc", bgc_id=bgc_id))

    return render_template(
        "edit/biosynth_class_specific.html",
        form=form,
        b_class=b_class,
        bgc_id=bgc_id,
        is_reviewer=current_user.has_role("reviewer"),
    )


@bp_edit.route("/<bgc_id>/biosynth/paths", methods=["GET", "POST"])
@login_required
def edit_biosynth_paths(bgc_id: str) -> str | response.Response:
    """Form to enter biosynthetic path information

    Args:
        bgc_id (str): BGC identifier

    Returns:
        str | Response: rendered template or redirect to edit_biosynth overview
    """
    if not is_valid_bgc_id(bgc_id):
        return abort(403, "Invalid existing entry!")

    if not request.form:
        form = FormCollection.paths(
            MultiDict(Storage.read_data(bgc_id).get("Biosynth_paths"))
        )
    else:
        form = FormCollection.paths(request.form)

    if request.method == "POST" and form.validate():
        # TODO: save to db
        Storage.save_data(bgc_id, "Biosynth_paths", request.form)
        flash("Submitted biosynthetic path information!")
        return redirect(url_for("edit.edit_biosynth", bgc_id=bgc_id))

    return render_template(
        "edit/biosynth_paths.html",
        bgc_id=bgc_id,
        form=form,
        is_reviewer=current_user.has_role("reviewer"),
    )


@bp_edit.route("/<bgc_id>/biosynth/modules", methods=["GET", "POST"])
@login_required
def edit_biosynth_modules(bgc_id: str) -> str | response.Response:
    """Form to enter biosynthetic module information

    Args:
        bgc_id (str): BGC identifier

    Returns:
        str | Response: rendered template or redirect to edit_biosynth overview
    """
    if not is_valid_bgc_id(bgc_id):
        return abort(403, "Invalid existing entry!")

    if not request.form:
        form = FormCollection.modules(
            MultiDict(Storage.read_data(bgc_id).get("Biosynth_modules"))
        )
    else:
        form = FormCollection.modules(request.form)

    if request.method == "POST" and form.validate():
        # TODO: save to db
        Storage.save_data(bgc_id, "Biosynth_modules", request.form)
        flash("Submitted biosynthetic module information!")
        return redirect(url_for("edit.edit_biosynth", bgc_id=bgc_id))

    return render_template(
        "edit/biosynth_modules.html",
        bgc_id=bgc_id,
        form=form,
        is_reviewer=current_user.has_role("reviewer"),
    )


@bp_edit.route("/<bgc_id>/tailoring", methods=["GET", "POST"])
@login_required
def edit_tailoring(bgc_id: str) -> str | response.Response:
    """Form to enter tailoring enzyme information

    Args:
        bgc_id (str): BGC identifier

    Returns:
        str | Response: rendered template or redirect to edit_bgc overview
    """
    if not is_valid_bgc_id(bgc_id):
        return abort(403, "Invalid existing entry!")

    if not request.form:
        form = FormCollection.tailoring(
            MultiDict(Storage.read_data(bgc_id).get("Tailoring"))
        )
    else:
        form = FormCollection.tailoring(request.form)

    if request.method == "POST" and form.validate():
        Storage.save_data(bgc_id, "Tailoring", request.form)
        flash("Submitted tailoring information!")
        return redirect(url_for("edit.edit_bgc", bgc_id=bgc_id))

    return render_template(
        "edit/tailoring.html",
        bgc_id=bgc_id,
        form=form,
        is_reviewer=current_user.has_role("reviewer"),
    )


@bp_edit.route("/<bgc_id>/annotation", methods=["GET", "POST"])
@login_required
def edit_annotation(bgc_id: str) -> str | response.Response:
    """Form to enter gene annotation information

    Args:
        bgc_id (str): BGC identifier

    Returns:
        str | Response: rendered template or redirect to edit_bgc overview
    """
    if not is_valid_bgc_id(bgc_id):
        return abort(403, "Invalid existing entry!")

    if not request.form:
        form = FormCollection.annotation(
            MultiDict(Storage.read_data(bgc_id).get("Annotation"))
        )
    else:
        form = FormCollection.annotation(request.form)

    if request.method == "POST" and form.validate():
        Storage.save_data(bgc_id, "Annotation", request.form)
        flash("Submitted annotation information!")
        return redirect(url_for("edit.edit_bgc", bgc_id=bgc_id))

    return render_template(
        "edit/annotation.html",
        bgc_id=bgc_id,
        form=form,
        is_reviewer=current_user.has_role("reviewer"),
    )


@bp_edit.route("/add_field", methods=["POST"])
def add_field() -> str:
    """Render an additional field as a subform

    Whenever a FieldList triggers this request for the addition of an entry to the list,
    use the id of the trigger to determine the where to append the entry.

    Examples:
        trigger id = 'structures'  -->   form.structures.append_entry()
        trigger id = 'enzymes-0-enzyme-0-auxiliary_enzymes'  -->
                        form.enzymes[0].enzyme[0].auxiliary_enzymes.append_entry()

    Returns:
        str: appended field rendered as a subform template string
    """
    # find directions to field to append an entry to
    directions = request.headers["Hx-Trigger"].split("-")
    # get origin of request to determine which form to use
    formname = Path(request.referrer).name
    curr = getattr(FormCollection, formname)(request.form)

    # sequentially traverse form fields
    i = 0
    while i + 2 < len(directions):
        subform, subform_idx = directions[i : i + 2]
        curr = getattr(curr, subform)[int(subform_idx)]
        i += 2

    # until we reach the final field that issued the request
    final = getattr(curr, directions[i])
    final.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        {{m.simple_divsubform(field, deletebtn=True)}}""",
        field=final[-1],
    )


@bp_edit.route("/get_references", methods=["POST"])
def get_references() -> str:
    bgc_id_match = re.search("edit/([^/]+)/", request.referrer)
    if bgc_id_match is not None:
        bgc_id = bgc_id_match.group(1)

    li = (
        lambda val: f"<li id={val} hx-post='/edit/append_reference' hx-target='previous input' hx-swap='outerHTML' hx-trigger='mousedown'>{val}</li>"
    )
    options = (
        "<span class='text-muted form-text'>Known references for this entry:</span>"
    )
    for ref in ReferenceUtils.collect_references(bgc_id):
        options += li(ref)
    return Markup(options)


@bp_edit.route("/append_reference", methods=["POST"])
def append_reference():

    target = request.headers.get("Hx-Target")
    current_val = request.values.get(target)
    new_ref = request.headers.get("Hx-Trigger")

    new_value = ReferenceUtils.append_ref(current_val, new_ref)
    return Markup(
        f"<input class='form-control' {html_params(value=new_value, id=target, name=target)}>"
    )
