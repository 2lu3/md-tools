"""step manager.

期待するプロジェクト構造

01, 02, 03, ... は step の番号。2桁とは限らない
(project_root)/
├── data/
│   ├── 01_(step name)/
│   ├── 02_(step name)/
│   ├── 03_(step name)/
│   └── ...
├── src/
│   ├── (package name)/
│   │   ├── 01_(step name)/
│   │   │   ├── template/
│   │   ├── 02_(step name)/
│   │   ├── 03_(step name)/
│   │   └── ...
├── template/
├── pyproject.toml
"""

import re
import shutil
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import toml
from alive_progress import alive_bar


@dataclass
class Step:
    """Represent one numbered project step."""

    idx: int
    name: str
    _project_root: Path
    _package_name: str

    @property
    def src_dir(self) -> Path:
        """Return the source directory for this step."""
        return (
            self._project_root
            / "src"
            / self._package_name
            / f"{self.idx:02d}_{self.name}"
        )

    @property
    def data_dir(self) -> Path:
        """Return the data directory for this step."""
        return self._project_root / "data" / f"{self.idx:02d}_{self.name}"

    @property
    def template_dir(self) -> Path:
        """Return the local template directory for this step."""
        return self.src_dir / "template"

    @property
    def global_template_dir(self) -> Path:
        """Return the package-level template directory for this step."""
        return self._project_root / "src" / self._package_name / "template"


class ProjectLayout:
    """Resolve project root, package name, and current or previous steps.

    steps はディレクトリ名の番号(例: 01_foo の 1)をキーとする dict。
    ・0 始まりを仮定しない。
    ・連番を前提にする
    """

    def __init__(self, start: Path | None = None) -> None:
        """Initialize project layout metadata from a start path."""
        self.project_root = self._find_project_root(start)
        self.package_name = self._find_package_name()
        self.steps: dict[int, Step] = self._list_steps()

    def current_step(self) -> Step:
        """Return the step containing the current working directory."""
        cwd = Path.cwd().resolve()
        for step in self.steps.values():
            if cwd.is_relative_to(step.src_dir.resolve()):
                return step
        msg = f"Current step not found from cwd: {cwd}"
        raise ValueError(msg)

    def previous_step(self) -> Step:
        """Return the step immediately before the current step."""
        previous_idx = self.current_step().idx - 1
        return self.steps[previous_idx]

    def get_step_by_index(self, index: int) -> Step:
        """Return a step by its numeric index."""
        return self.steps[index]

    @property
    def global_template_dir(self) -> Path:
        """Return the project-level template directory."""
        return self.project_root / "src" / "template"

    def _find_package_name(self) -> str:
        # pyproject.toml の [project] セクションの name を取得
        with (self.project_root / "pyproject.toml").open() as f:
            toml_data = toml.load(f)
        return toml_data["project"]["name"]

    def _list_steps(self) -> dict[int, Step]:
        """List step directories matching src/package_name/(number)_(step_name)/."""
        steps: dict[int, Step] = {}
        steps_root = self.project_root / "src" / self.package_name

        for step_dir in steps_root.iterdir():
            if not step_dir.is_dir():
                continue

            match = re.match(r"^(\d+)_(.+)$", step_dir.name)
            if match is None:
                continue

            index = int(match.group(1))
            name = match.group(2)
            if steps.get(index) is not None:
                msg = f"Step index {index} already exists"
                raise ValueError(msg)
            steps[index] = Step(
                idx=index,
                name=name,
                _project_root=self.project_root,
                _package_name=self.package_name,
            )

        if not steps:
            msg = f"Step directories not found in {steps_root}"
            raise FileNotFoundError(msg)

        return dict(sorted(steps.items()))

    def _find_project_root(self, start: Path | None = None) -> Path:
        """Return the first parent directory containing pyproject.toml."""
        origin = (start or Path.cwd()).resolve()
        p = origin
        while True:
            if (p / "pyproject.toml").is_file():
                return p
            if p.parent == p:
                msg = f"pyproject.toml not found from {origin}"
                raise FileNotFoundError(msg)
            p = p.parent


class StepManager:
    """Manage data directories and file handoff between project steps."""

    def __init__(
        self,
        condition_ids: Sequence[str] = (),
        layout: ProjectLayout | None = None,
    ) -> None:
        """Initialize the manager with condition IDs and an optional layout."""
        self._condition_ids: tuple[str, ...] = tuple(condition_ids)
        self._layout: ProjectLayout = layout or ProjectLayout()

    def init_data_dir(
        self,
        subdirs: Sequence[str] = ("out", "in", "script", "next_in"),
        non_ignore_files: Sequence[str] = ("*.dvc", ".gitignore"),
    ) -> None:
        """Create a fresh data directory for the current step."""
        current_step = self._layout.current_step()
        shutil.rmtree(current_step.data_dir, ignore_errors=True)
        current_step.data_dir.mkdir(parents=True)

        for subdir in subdirs:
            (current_step.data_dir / subdir).mkdir()

        with (current_step.data_dir / ".gitignore").open("w") as f:
            f.write("*\n")
            f.writelines(
                f"!{non_ignore_file}\n" for non_ignore_file in non_ignore_files
            )
        with (current_step.data_dir / ".gitkeep").open("w") as f:
            f.write("")

    def copy_files(self) -> None:
        """Copy previous step next_in files to the current step input directory."""
        src_dir = self._layout.previous_step().data_dir / "next_in"
        dst_dir = self._layout.current_step().data_dir / "in"

        shutil.rmtree(dst_dir, ignore_errors=True)
        shutil.copytree(src_dir, dst_dir)

    def prepare_next_in(self, copy_files: list[tuple[str, str, list[str]]]) -> None:
        """Prepare current step output files as next step input files."""
        with alive_bar(
            len(self._condition_ids) * len(copy_files), title="Copying files"
        ) as bar:
            for src_dir, suffix, extensions in copy_files:
                for condition_id in self._condition_ids:
                    for extension in extensions:
                        if suffix != "":
                            src_file = (
                                self.cur_step.data_dir
                                / src_dir
                                / f"{condition_id}_{suffix}.{extension}"
                            )
                        else:
                            src_file = (
                                self.cur_step.data_dir
                                / src_dir
                                / f"{condition_id}.{extension}"
                            )
                        shutil.copyfile(
                            src_file,
                            self.cur_step.data_dir
                            / "next_in"
                            / f"{condition_id}.{extension}",
                        )
                    bar()

    @property
    def cur_step(self) -> Step:
        """Return the current step."""
        return self._layout.current_step()

    @property
    def prev_step(self) -> Step:
        """Return the previous step."""
        return self._layout.previous_step()

    def step(self, index: int) -> Step:
        """Return a step by index."""
        return self._layout.get_step_by_index(index)
