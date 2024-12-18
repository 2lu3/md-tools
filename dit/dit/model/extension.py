from dit.model.config import Configuration
from .git import Git

class Extension:
    """dvcで管理されるべき拡張子を示す"""
    def __init__(self):
       self.extensions: list[str] = self._load()

    def add(self, pattern: str):
        self.extensions.append(pattern)
        self._save()

    def remove(self, pattern: str):
        assert pattern in self.extensions
        self.extensions.remove(pattern)
        self._save()


    def is_match(self, path: str):
        for ext in self.extensions:
            if path.endswith(ext):
                return True
        return False


    def _load(self):
        git = Git()
        config = Configuration(git.root_dir())
        extensions: list[str] = config.load_config().get("extensions", [])
        assert type(extensions) == list
        return extensions

    def _save(self):
        git = Git()
        config = Configuration(git.root_dir())
        config.update({"extensions": self.extensions})
