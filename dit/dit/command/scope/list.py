import argparse
from dit.model.scope import Scope
from natsort import natsorted

def list_scope():
    scope = Scope()

    directories = scope.directories

    for directory in natsorted(directories):
        print(directory)



def to_command():
    parser = argparse.ArgumentParser()

    args = parser.parse_args()

    list_scope()

