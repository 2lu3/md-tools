"""Render Jinja templates for configuration files."""

from collections.abc import Mapping
from pathlib import Path

from jinja2 import Environment, FileSystemLoader


class TemplateRenderer:
    """Render templates from a configured template directory."""

    def __init__(
        self,
        params: Mapping[str, object],
        template_dir: str = "template",
    ) -> None:
        """Initialize the template renderer.

        Args:
            params: Template parameters.
            template_dir: Directory containing templates.
        """
        self._params = params
        self._template_dir = template_dir
        self._env = Environment(
            loader=FileSystemLoader(self._template_dir),
            autoescape=False,  # noqa: S701 - templates render non-HTML inputs.
        )

    def render(self, template_filename: str, output_path: str) -> None:
        """Render a template to a file.

        Args:
            template_filename: Template file name.
            output_path: Output file path.
        """
        template = self._env.get_template(template_filename)

        rendered = template.render(**self._params)

        with Path(output_path).open("w") as f:
            f.write(rendered)


def render2file(
    output_path: str,
    template_path: str,
    params: Mapping[str, object],
) -> None:
    """Render a template to a file.

    Args:
        output_path: Output file path.
        template_path: Template file name.
        params: Template parameters.
    """
    renderer = TemplateRenderer(params)
    renderer.render(template_path, output_path)
