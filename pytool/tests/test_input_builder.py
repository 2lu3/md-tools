import os
import shutil

from pytool.input_builder import write_template

def test_input_builder():

    shutil.rmtree('sandbox', ignore_errors=True)
    os.makedirs('sandbox')
    os.chdir('sandbox')

    write_template('min')

    os.chdir('..')
    shutil.rmtree('sandbox', ignore_errors=True)
