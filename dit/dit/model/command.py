from abc import ABC, abstractmethod


class AbstCommand(ABC):
    @abstractmethod
    def register_subparser(self, subparser):
        pass

    @abstractmethod
    def handler(self, args):
        pass

