from typing import Generator, Optional
from abc import ABCMeta, abstractmethod
import os

from pytool.config.render_template import TemplateRenderer


class ProjectBuilder(metaclass=ABCMeta):
    """ProjectBuilder."""

    def __init__(
        self, project_name: str, is_benchmark: bool, is_clean_dirs: bool = False
    ):
        self.project_name = project_name
        self.is_benchmark = is_benchmark
        self.is_clean_dirs = is_clean_dirs

    @abstractmethod
    def project_params(self) -> Generator[dict, None, None]:
        pass

    @abstractmethod
    def benchmark_params(self) -> Generator[dict, None, None]:
        pass

    @abstractmethod
    def create_project(self, param: dict):
        pass

    def create_benchmark(self, param: dict):
        self.create_project(param)

    def write_input(
        self,
        project_dir: str,
        mode: str,
        index: int,
        template_file: Optional[str] = None,
        **kwargs,
    ):
        """write_input.

        Args:
            project_dir (str): project_dir
            mode (str): mode
            index (int): index
            template_file (Optional[str]): template_file
            kwargs:
        """
        params = self._local2params(locals(), kwargs)

        if template_file is None:
            template_file = f"{mode}{index}.inp"

        renderer = TemplateRenderer(params)
        renderer.render(
            template_file, os.path.join(project_dir, "inp", f"{mode}{index}.inp")
        )

    def write_job(
        self,
        project_dir: str,
        mode: str,
        index: int,
        template_file: str = "job.sh",
        **kwargs,
    ):
        """write_job.

        Args:
            project_dir (str): project_dir
            mode (str): mode
            index (int): index
            template_file (str): template_file
            kwargs:
        """
        params = self._local2params(locals(), kwargs)

        renderer = TemplateRenderer(params)
        renderer.render(
            template_file, os.path.join(project_dir, f"job_{mode}{index}.sh")
        )
        os.chmod(os.path.join(project_dir, f"job_{mode}{index}.sh"), 0o755)

    def write_submit(self, project_dir: str, **kwargs):
        """write_submit.

        Args:
            project_dir (str): project_dir
            kwargs:
        """
        params = self._local2params(locals(), kwargs)

        renderer = TemplateRenderer(params)
        renderer.render("submit.sh", os.path.join(project_dir, "submit.sh"))
        os.chmod(os.path.join(project_dir, "submit.sh"), 0o755)

    @staticmethod
    def write_submit_projects(project_dirs: list[str]):
        TemplateRenderer({"projects": project_dirs}).render(
            "submit_projects.sh", "submit_projects.sh"
        )
        os.chmod("submit_projects.sh", 0o755)

    @staticmethod
    def write_submit_benchmarks(is_benchmark: bool):
        TemplateRenderer({}).render("submit_benchmarks.sh", "submit_benchmarks.sh")
        os.chmod("submit_benchmarks.sh", 0o755)

    def build(self):
        if self.is_benchmark:
            os.makedirs("benchmarks", exist_ok=True)
            for param in self.benchmark_params():
                self.create_benchmark(param)
        else:
            for param in self.project_params():
                self.create_project(param)

    def _local2params(self, params: dict, kwargs: dict):
        """locals()で取得したローカル変数とkwargsを合成する
        ローカル変数のうちselfとkwargsは除外する

        Args:
            params (dict): params
            kwargs (dict): kwargs
        """
        del params["self"]
        del params["kwargs"]
        params.update(kwargs)

        return params
