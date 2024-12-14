import os

class Git:
    def __init__(self):
        pass

    def root_dir(self):
        current_dir = os.getcwd()
        while not os.path.exists(os.path.join(current_dir, '.git')):
            current_dir = os.path.dirname(current_dir)
        return current_dir
