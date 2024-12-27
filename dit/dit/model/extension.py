from dit.model.config import Configuration
from .git import Git

class Extension:
    extensions: set[str] = set()

    """dvcで管理されるべき拡張子を示す"""
    def __init__(self):
        self._load()

    def add(self, pattern: str):
        Extension.extensions.add(pattern)
        self._save()

    def remove(self, pattern: str):
        assert pattern in Extension.extensions
        Extension.extensions.remove(pattern)
        self._save()


    def is_match(self, path: str):
        for ext in Extension.extensions:
            if path.endswith(ext):
                return True
        return False


    def _load(self):
        git = Git()
        config = Configuration()
        extensions: list[str] = config.load_config().get("extensions", [])
        assert type(extensions) == list
        Extension.extensions = set(extensions)

    def _save(self):
        git = Git()
        config = Configuration()
        config.update({"extensions": list(Extension.extensions)})
