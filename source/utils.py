from pathlib import Path
import sys
import numpy as np

def get_path_root():
    if sys.platform == 'linux':
        path_root = Path(sys.executable).parent
    elif sys.platform == 'win32':
        path_root = Path(sys.executable).parent
    else:
        raise ValueError('Unsupported OS.')
    
    return path_root


# Helper functions
def locate_file(from_path : Path = None, query_path = None, query_name = "query file"):

    if not from_path or not query_path:
        raise ValueError("Must specify from_path and query_path")


    possible_path = list(from_path.rglob(query_path))

    if not possible_path:
        raise FileNotFoundError(f"Cannot find {query_name} from {from_path} by {query_path}")

    return_which = (
        f"using {query_name} at:\n"
        f"{possible_path[0]}\n"
    )
    print(return_which)

    return possible_path[0]

def l2_norm(x):
    return np.sqrt(np.sum(np.asarray(x) ** 2))
