import PySimpleGUI as sg
from dit.model.scope import Scope

class GUI:
    def __init__(self):
        self.window = sg.Window("dit", self._layout())

        while True:
            event, values = self.window.read()
            if event == sg.WIN_CLOSED:
                break
            if self.loop(event, values):
                break
        self.window.close()

    def loop(self, event, values):

        return False

    def _layout(self):
        return [
                []

                ]
