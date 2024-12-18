import os
import json
from loguru import logger

class Configuration:
    def __init__(self, project_name: str):
        self.project_name = project_name
        if not os.path.exists(self._root_dir()):
            os.makedirs(self._root_dir())

        if not os.path.exists(self._config_file_path()):
            self.save({})

    def _root_dir(self):
        userhome_dir = os.path.expanduser('~')
        return os.path.join(userhome_dir, '.config', 'dit')

    def load_config(self):
        try:
            with open(self._config_file_path(), 'r') as f:
                return json.loads(f.read())
        except FileNotFoundError:
            logger.error(f'Config file not found: {self._config_file_path()}')
            return {}
        except json.JSONDecodeError:
            logger.error(f'Config file is not a valid JSON file: {self._config_file_path()}')
            return {}

    def save(self, config: dict):
        with open(self._config_file_path(), 'w') as f:
            f.write(json.dumps(config))

    def update(self, config: dict):
        current_config = self.load_config()
        current_config.update(config)
        self.save(current_config)

    def _config_file_path(self) -> str:
        return os.path.join(self._root_dir(), f'{self.project_name}.json')

