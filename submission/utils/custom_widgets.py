""" Collection of custom widget classes used throughout the submission system """

from typing import Optional

from flask import url_for
from wtforms import Field, SelectFieldBase, widgets
from markupsafe import Markup


class StringFieldAddBtn(widgets.TextInput):
    def __init__(self, label="add", render_kw={}) -> None:
        super(StringFieldAddBtn, self).__init__()
        self.label = label
        self.render_kw = render_kw

    def __call__(self, field, **kwargs):
        orig = super(StringFieldAddBtn, self).__call__(field, **kwargs)
        add_btn = f"<button class='btn btn-light' {self.html_params(**self.render_kw)}>{self.label}</button>"
        return orig + Markup(add_btn)


class FieldListAddBtn(widgets.SubmitInput):
    def __init__(
        self,
        input_type: Optional[str] = None,
        label="",
        render_kw={
            "formnovalidate": True,
            "hx-post": "/edit/add_field",
            "hx-swap": "beforebegin",
        },
    ) -> None:
        super().__init__(input_type)
        self.label = label
        self.render_kw = render_kw

    def __call__(self, field: Field, **kwargs: object) -> Markup:
        kwargs.setdefault("id", field.id)
        kwargs.setdefault("type", self.input_type)
        return Markup(
            f"<button class='btn btn-light' {self.html_params(**kwargs, **self.render_kw)}>{self.label}</button>"
        )


class MultiTextInput(widgets.TextInput):
    def __init__(
        self, input_type: Optional[str] = None, number: int = 1, render_kw: list[dict] = []
    ) -> None:
        super().__init__(input_type)
        self.number = number
        self.render_kw = render_kw

    def __call__(self, field: Field, **kwargs: object) -> Markup:
        kwargs.setdefault("id", field.id)
        kwargs.setdefault("type", self.input_type)

        inputs = Markup("<div>")
        for field_idx in range(self.number):
            render_kw = self.render_kw[field_idx]
            inputs += Markup("<input %s>" % self.html_params(**render_kw, **kwargs))
        inputs += Markup("</div>")
        return inputs


class SelectDefault(widgets.Select):
    def __init__(
        self,
        multiple: bool = False,
        default: Optional[tuple[str, str]] = ("", " --- Select --- "),
        **kwargs,
    ) -> None:
        super().__init__(multiple)
        self.default = default
        self.kwargs = kwargs

    def __call__(self, field: SelectFieldBase, **kwargs: object) -> Markup:
        select = super().__call__(field, **kwargs)
        if self.default:
            val, label = self.default
            insert_idx = select.find(">") + 1
            default_option = Markup(
                f"<option selected {widgets.html_params(value=val, **kwargs, **self.kwargs)}>{label}</option>"
            )
            select = select[:insert_idx] + default_option + select[insert_idx:]
        return select


class TextInputIndicator(widgets.TextInput):
    def __init__(self, input_type: Optional[str] = None) -> None:
        super().__init__(input_type)

    def __call__(self, field: Field, **kwargs: object) -> Markup:
        spinner = Markup(
            f'<img id="spinner" class="htmx-indicator" style="display:inline-block;margin:5px" src="{url_for("static", filename="img/wifi-fade.svg")}" />'
        )
        kwargs.setdefault("style", "width:95%;display:inline-block")
        return super().__call__(field, **kwargs) + spinner


class StructureInput(widgets.TextInput):
    def __init__(
        self,
        input_type: Optional[str] = None,
        render_kw={
            "hx-post": "/edit/render_smiles",
            "hx-trigger": "change, load",
            "hx-swap": "innerHTML",
            "hx-target": "next .struct",
        },
    ) -> None:
        super().__init__(input_type)
        self.render_kw = render_kw

    def __call__(self, field: Field, **kwargs: object) -> Markup:

        return super().__call__(field, **kwargs, **self.render_kw) + Markup(
            "<div class='struct'></div>"
        )


class TextInputWithSuggestions(widgets.TextInput):
    def __init__(
        self, input_type: Optional[str] = None, post_url: Optional[str] = None
    ) -> None:
        super().__init__(input_type)
        self.render_kw = {
            "hx-post": post_url,
            "hx-trigger": "load",
            "hx-target": "next ul",
            "hx-swap": "innerHTML",
        }

    def __call__(self, field: Field, **kwargs: object) -> Markup:
        obj = Markup("<div class='suggestions'>")
        obj += super().__call__(field, **kwargs, **self.render_kw)
        obj += Markup("<ul class='form-control'></ul></div>")
        return obj
