from jinja2 import Environment, FileSystemLoader
import os


class TemplateRenderer:
    def __init__(self, params: dict, template_dir: str = "template"):
        self._params = params
        self._template_dir = template_dir
        self._env = Environment(loader=FileSystemLoader(self._template_dir))

    def render(self, template_filename: str, output_path: str):
        template = self._env.get_template(template_filename)

        rendered = template.render(**self._params)

        with open(output_path, "w") as f:
            f.write(rendered)

def render2file(output_path: str, template_path: str, params: dict):
    renderer = TemplateRenderer(params)
    renderer.render(template_path, output_path)
