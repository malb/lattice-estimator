import argparse
import os
import subprocess

from nbconvert import HTMLExporter
from nbconvert.preprocessors import ClearOutputPreprocessor, ExecutePreprocessor
from nbconvert.preprocessors import CellExecutionError
import nbformat

def process_notebook(notebook_filename, html_directory = 'notebook-html', execute=True):
    '''Checks if an IPython notebook runs without error from start to finish. If so, writes the notebook to HTML (with outputs) and overwrites the .ipynb file (without outputs).
    '''
    with open(notebook_filename) as f:
        nb = nbformat.read(f, as_version=4)
        
    clear = ClearOutputPreprocessor()
    
    try:
        if execute:
            # Check that the notebook runs
            ep = ExecutePreprocessor(timeout=600, kernel_name='python3')
            ep.preprocess(nb, {'metadata': {'path': ''}})

        msg = ''

    except CellExecutionError:
        out = None
        msg = f'\n  Error executing the notebook "{notebook_filename}".\n'
        msg += f'  See notebook "{notebook_filename}" for the traceback.'
        
    # Clear notebook outputs and save as .ipynb
    cleared = clear.preprocess(nb, {})
    with open(notebook_filename, mode='w', encoding='utf-8') as f:
        nbformat.write(nb, f)
         
    print(f"Processed {notebook_filename}{msg}")

    return


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='read some notebok files')
    parser.add_argument('notebooks',
                        metavar='Notebooks', 
                        type=str, 
                        nargs='+',
                        help='notebooks')

    args = parser.parse_args()

    notebooks = args.notebooks

    for fn in notebooks:
        if not fn.endswith('.ipynb'):
            print(f'Error: file {fn} is not an IPython notebook.')
            raise
        
    for fn in notebooks:
        process_notebook(fn, execute=False)
    