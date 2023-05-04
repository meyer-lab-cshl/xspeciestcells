### Adapted from https://github.com/ventolab/CellphoneDB/blob/master/notebooks/T00_DownloadDB.ipynb

import pandas as pd
import sys
import os

from IPython.display import HTML, display
from cellphonedb.utils import db_releases_utils
from cellphonedb.utils import db_utils

pd.set_option('display.max_columns', 100)

# Define our base directory for the analysis
directory = '/Users/hmeyer/projects/' + \
            '20220809_Thymic-iNKT-CrossSpecies/cellphonedb'

# Version of the databse
cpdb_version = 'v4.1.0'

# Path where the input files to generate the database are located
cpdb_target_dir = os.path.join(directory, "db", cpdb_version)

# Download database
db_utils.download_database(cpdb_target_dir, cpdb_version)
