{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Liver Model Construction: Notebook-Glycolysis_Gluconeogenesis-trial"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Setup workflow"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Import packages"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "MASSpy version: 0.1.1\n"
     ]
    }
   ],
   "source": [
    "import os\n",
    "import warnings\n",
    "from cobra.io.json import load_json_model as load_json_cobra_model\n",
    "import escher\n",
    "import mass\n",
    "import numpy as np\n",
    "import pandas as pd\n",
    "import sympy as sy\n",
    "from cobra import Model, Reaction, Metabolite\n",
    "import cobra.test\n",
    "from os.path import join\n",
    "from mass.util import qcqa_model\n",
    "from cobra import DictList\n",
    "from mass import (\n",
    "    MassConfiguration, MassMetabolite, MassModel,\n",
    "    MassReaction, Simulation, UnitDefinition)\n",
    "from mass.io.json import save_json_model as save_json_mass_model\n",
    "from mass.visualization import plot_comparison, plot_time_profile\n",
    "mass_config = MassConfiguration()\n",
    "from six import iteritems\n",
    "\n",
    "mass_config.irreversible_Keq = float(\"inf\")\n",
    "import matplotlib.pyplot as plt\n",
    "\n",
    "print(\"MASSpy version: {0}\".format(mass.__version__))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [],
   "source": [
    "# This takes in the file name (filepath) and sheet name to then extract data to then \n",
    "def load_data(filepath, sheet_name):\n",
    "    \"\"\"Load Liver data from an excel sheet\"\"\"\n",
    "    df = pd.read_excel(\n",
    "        io=filepath,\n",
    "        sheet_name=sheet_name,\n",
    "        index_col=0)\n",
    "    df = df.drop(\"Additional Notes\", axis=1)\n",
    "    df = df.drop(\"Name\", axis=1)\n",
    "    return df"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {},
   "outputs": [
    {
     "ename": "XLRDError",
     "evalue": "Excel xlsx file; not supported",
     "output_type": "error",
     "traceback": [
      "\u001b[0;31m---------------------------------------------------------------------------\u001b[0m",
      "\u001b[0;31mXLRDError\u001b[0m                                 Traceback (most recent call last)",
      "\u001b[0;32m<ipython-input-3-ef4d4b1a75f0>\u001b[0m in \u001b[0;36m<module>\u001b[0;34m\u001b[0m\n\u001b[0;32m----> 1\u001b[0;31m \u001b[0mconc_df\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mload_data\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m\"1.2-ma-data-collection.xlsx\"\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;34m\"Concentrations\"\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m      2\u001b[0m \u001b[0mconc_df\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m<ipython-input-2-051e3e8201a9>\u001b[0m in \u001b[0;36mload_data\u001b[0;34m(filepath, sheet_name)\u001b[0m\n\u001b[1;32m      5\u001b[0m         \u001b[0mio\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0mfilepath\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m      6\u001b[0m         \u001b[0msheet_name\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0msheet_name\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m----> 7\u001b[0;31m         index_col=0)\n\u001b[0m\u001b[1;32m      8\u001b[0m     \u001b[0mdf\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mdf\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mdrop\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m\"Additional Notes\"\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0maxis\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0;36m1\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m      9\u001b[0m     \u001b[0mdf\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mdf\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mdrop\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m\"Name\"\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0maxis\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0;36m1\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m/usr/local/lib/python3.7/site-packages/pandas/util/_decorators.py\u001b[0m in \u001b[0;36mwrapper\u001b[0;34m(*args, **kwargs)\u001b[0m\n\u001b[1;32m    294\u001b[0m                 )\n\u001b[1;32m    295\u001b[0m                 \u001b[0mwarnings\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mwarn\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mmsg\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mFutureWarning\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mstacklevel\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0mstacklevel\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 296\u001b[0;31m             \u001b[0;32mreturn\u001b[0m \u001b[0mfunc\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m*\u001b[0m\u001b[0margs\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;34m**\u001b[0m\u001b[0mkwargs\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    297\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    298\u001b[0m         \u001b[0;32mreturn\u001b[0m \u001b[0mwrapper\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m/usr/local/lib/python3.7/site-packages/pandas/io/excel/_base.py\u001b[0m in \u001b[0;36mread_excel\u001b[0;34m(io, sheet_name, header, names, index_col, usecols, squeeze, dtype, engine, converters, true_values, false_values, skiprows, nrows, na_values, keep_default_na, na_filter, verbose, parse_dates, date_parser, thousands, comment, skipfooter, convert_float, mangle_dupe_cols)\u001b[0m\n\u001b[1;32m    302\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    303\u001b[0m     \u001b[0;32mif\u001b[0m \u001b[0;32mnot\u001b[0m \u001b[0misinstance\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mio\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mExcelFile\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 304\u001b[0;31m         \u001b[0mio\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mExcelFile\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mio\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mengine\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0mengine\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    305\u001b[0m     \u001b[0;32melif\u001b[0m \u001b[0mengine\u001b[0m \u001b[0;32mand\u001b[0m \u001b[0mengine\u001b[0m \u001b[0;34m!=\u001b[0m \u001b[0mio\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mengine\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    306\u001b[0m         raise ValueError(\n",
      "\u001b[0;32m/usr/local/lib/python3.7/site-packages/pandas/io/excel/_base.py\u001b[0m in \u001b[0;36m__init__\u001b[0;34m(self, path_or_buffer, engine)\u001b[0m\n\u001b[1;32m    865\u001b[0m         \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0m_io\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mstringify_path\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mpath_or_buffer\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    866\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 867\u001b[0;31m         \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0m_reader\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0m_engines\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0mengine\u001b[0m\u001b[0;34m]\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0m_io\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    868\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    869\u001b[0m     \u001b[0;32mdef\u001b[0m \u001b[0m__fspath__\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mself\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m/usr/local/lib/python3.7/site-packages/pandas/io/excel/_xlrd.py\u001b[0m in \u001b[0;36m__init__\u001b[0;34m(self, filepath_or_buffer)\u001b[0m\n\u001b[1;32m     20\u001b[0m         \u001b[0merr_msg\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;34m\"Install xlrd >= 1.0.0 for Excel support\"\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m     21\u001b[0m         \u001b[0mimport_optional_dependency\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m\"xlrd\"\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mextra\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0merr_msg\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m---> 22\u001b[0;31m         \u001b[0msuper\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0m__init__\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mfilepath_or_buffer\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m     23\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m     24\u001b[0m     \u001b[0;34m@\u001b[0m\u001b[0mproperty\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m/usr/local/lib/python3.7/site-packages/pandas/io/excel/_base.py\u001b[0m in \u001b[0;36m__init__\u001b[0;34m(self, filepath_or_buffer)\u001b[0m\n\u001b[1;32m    351\u001b[0m             \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mbook\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mload_workbook\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mfilepath_or_buffer\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    352\u001b[0m         \u001b[0;32melif\u001b[0m \u001b[0misinstance\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mfilepath_or_buffer\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mstr\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 353\u001b[0;31m             \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mbook\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mload_workbook\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mfilepath_or_buffer\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    354\u001b[0m         \u001b[0;32melif\u001b[0m \u001b[0misinstance\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mfilepath_or_buffer\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mbytes\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    355\u001b[0m             \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mbook\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mload_workbook\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mBytesIO\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mfilepath_or_buffer\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m/usr/local/lib/python3.7/site-packages/pandas/io/excel/_xlrd.py\u001b[0m in \u001b[0;36mload_workbook\u001b[0;34m(self, filepath_or_buffer)\u001b[0m\n\u001b[1;32m     35\u001b[0m             \u001b[0;32mreturn\u001b[0m \u001b[0mopen_workbook\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mfile_contents\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0mdata\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m     36\u001b[0m         \u001b[0;32melse\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m---> 37\u001b[0;31m             \u001b[0;32mreturn\u001b[0m \u001b[0mopen_workbook\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mfilepath_or_buffer\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m     38\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m     39\u001b[0m     \u001b[0;34m@\u001b[0m\u001b[0mproperty\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.7/site-packages/xlrd/__init__.py\u001b[0m in \u001b[0;36mopen_workbook\u001b[0;34m(filename, logfile, verbosity, use_mmap, file_contents, encoding_override, formatting_info, on_demand, ragged_rows, ignore_workbook_corruption)\u001b[0m\n\u001b[1;32m    168\u001b[0m     \u001b[0;31m# files that xlrd can parse don't start with the expected signature.\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    169\u001b[0m     \u001b[0;32mif\u001b[0m \u001b[0mfile_format\u001b[0m \u001b[0;32mand\u001b[0m \u001b[0mfile_format\u001b[0m \u001b[0;34m!=\u001b[0m \u001b[0;34m'xls'\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 170\u001b[0;31m         \u001b[0;32mraise\u001b[0m \u001b[0mXLRDError\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mFILE_FORMAT_DESCRIPTIONS\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0mfile_format\u001b[0m\u001b[0;34m]\u001b[0m\u001b[0;34m+\u001b[0m\u001b[0;34m'; not supported'\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    171\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    172\u001b[0m     bk = open_workbook_xls(\n",
      "\u001b[0;31mXLRDError\u001b[0m: Excel xlsx file; not supported"
     ]
    }
   ],
   "source": [
    "conc_df = load_data(\"1.2-ma-data-collection.xlsx\", \"Concentrations\")\n",
    "conc_df"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Defaulting to user installation because normal site-packages is not writeable\n",
      "Requirement already satisfied: cplex in ./.local/lib/python3.7/site-packages (20.1.0.1)\n",
      "\u001b[33mWARNING: You are using pip version 20.2.4; however, version 21.0.1 is available.\n",
      "You should consider upgrading via the '/usr/local/bin/python -m pip install --upgrade pip' command.\u001b[0m\n",
      "Note: you may need to restart the kernel to use updated packages.\n"
     ]
    }
   ],
   "source": [
    "pip install cplex"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Defaulting to user installation because normal site-packages is not writeable\n",
      "Requirement already satisfied: xlrd in ./.local/lib/python3.7/site-packages (2.0.1)\n",
      "\u001b[33mWARNING: You are using pip version 20.2.4; however, version 21.0.1 is available.\n",
      "You should consider upgrading via the '/usr/local/bin/python -m pip install --upgrade pip' command.\u001b[0m\n",
      "Note: you may need to restart the kernel to use updated packages.\n"
     ]
    }
   ],
   "source": [
    "pip install xlrd"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# from data_input_functions import( load_data)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Import functions "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "## load functions"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#### Set Configuration bounds"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "metadata": {},
   "outputs": [],
   "source": [
    "seed = int(4) ## starting point for randomised models, to ensure same output\n",
    "n_models = 100\n",
    "mass_config = MassConfiguration()\n",
    "mass_config.solver = \"cplex\"\n",
    "# mass_config.solver = \"glpk\"\n",
    "mass_config.bounds=(-1000,1000)\n",
    "## need to install cplex into the repository requirements"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#### Directory paths"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "metadata": {},
   "outputs": [],
   "source": [
    "model_dir = os.path.abspath(\"../mass_user/models\")\n",
    "maps_dir = os.path.abspath(\"../mass_user/maps\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "metadata": {},
   "outputs": [],
   "source": [
    "## Potential Conversion factors\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#### Add Enzyme Module data"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "metadata": {},
   "outputs": [],
   "source": [
    "## if adding enzyme modules, then add paths to include enzyme module data\n",
    "## add datasheet for data that may be needed"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Allow Escher to close without pop-up\n",
    "escher.rc['never_ask_before_quit'] = True"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Create Model of Glycolytic Network"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Load COBRA model"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "\n",
       "        <table>\n",
       "            <tr>\n",
       "                <td><strong>Name</strong></td>\n",
       "                <td>None</td>\n",
       "            </tr><tr>\n",
       "                <td><strong>Memory address</strong></td>\n",
       "                <td>0x07f0916daf690</td>\n",
       "            </tr><tr>\n",
       "                <td><strong>Number of metabolites</strong></td>\n",
       "                <td>178</td>\n",
       "            </tr><tr>\n",
       "                <td><strong>Number of reactions</strong></td>\n",
       "                <td>187</td>\n",
       "            </tr><tr>\n",
       "                <td><strong>Number of groups</strong></td>\n",
       "                <td>0</td>\n",
       "            </tr><tr>\n",
       "                <td><strong>Objective expression</strong></td>\n",
       "                <td>1.0*ATPM - 1.0*ATPM_reverse_5b752</td>\n",
       "            </tr><tr>\n",
       "                <td><strong>Compartments</strong></td>\n",
       "                <td>m, i, c, r, </td>\n",
       "            </tr>\n",
       "          </table>"
      ],
      "text/plain": [
       "<Model None at 0x7f0916daf690>"
      ]
     },
     "execution_count": 11,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "core_model=load_json_cobra_model(filename=os.path.join(model_dir,\"CoreModel.json\"))\n",
    "core_model"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Create MASS model from COBRA model"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Create MASS model"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "metadata": {},
   "outputs": [],
   "source": [
    "glycolysis = MassModel(\"Glycolysis\", \n",
    "                       array_type='DataFrame', \n",
    "                       dtype=np.int64)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### View  network to be built"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 13,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "application/vnd.jupyter.widget-view+json": {
       "model_id": "b51cf040d20841cfad502414babc5ce9",
       "version_major": 2,
       "version_minor": 0
      },
      "text/plain": [
       "Builder(highlight_missing=True, never_ask_before_quit=True)"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "### making sure core_model has all the reactions to start\n",
    "escher_builder = escher.Builder(\n",
    "    model=core_model,\n",
    "    map_json=os.path.join(\n",
    "        maps_dir, \".\".join((\n",
    "            glycolysis.id, \"map\", \"json\"))\n",
    "    ),\n",
    "    highlight_missing=True)\n",
    "\n",
    "escher_builder"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 14,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<strong><em>Optimal</em> solution with objective value 25.200</strong><br><div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>fluxes</th>\n",
       "      <th>reduced_costs</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>CSm</th>\n",
       "      <td>2.0</td>\n",
       "      <td>0.0</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>ACONTm</th>\n",
       "      <td>2.0</td>\n",
       "      <td>0.0</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>ICDHxm</th>\n",
       "      <td>2.0</td>\n",
       "      <td>0.0</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>AKGDm</th>\n",
       "      <td>2.0</td>\n",
       "      <td>0.0</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>SUCOASm</th>\n",
       "      <td>-2.0</td>\n",
       "      <td>-0.0</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>...</th>\n",
       "      <td>...</td>\n",
       "      <td>...</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>EX_nh4_c</th>\n",
       "      <td>0.0</td>\n",
       "      <td>4.0</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>EX_so3_c</th>\n",
       "      <td>-0.0</td>\n",
       "      <td>0.0</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>EX_etoh_c</th>\n",
       "      <td>-0.0</td>\n",
       "      <td>0.0</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>EX_glyc_3octa_c</th>\n",
       "      <td>-0.0</td>\n",
       "      <td>0.0</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>EX_fru_c</th>\n",
       "      <td>0.0</td>\n",
       "      <td>-50.4</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "<p>187 rows × 2 columns</p>\n",
       "</div>"
      ],
      "text/plain": [
       "<Solution 25.200 at 0x7f095954d490>"
      ]
     },
     "execution_count": 14,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "flux_solution = core_model.optimize()\n",
    "flux_solution\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "metadata": {},
   "outputs": [],
   "source": [
    "\n",
    "reaction = ['GLCter',\n",
    "           'G6PPer',\n",
    "           'G6Pter',\n",
    "           'GLPASE1',\n",
    "           'PGMT',\n",
    "            'EX_glygn2_c',\n",
    "           'GAPD',\n",
    "           'PGK',\n",
    "           'PGM',\n",
    "            'ENO',\n",
    "            'PEPCKm',\n",
    "            'PCm',\n",
    "            'LDH_L',\n",
    "           'ME2',\n",
    "            'PYK',\n",
    "           'CSm',\n",
    "           'SUCOASm',\n",
    "            'SUCD1m',\n",
    "            'FUMm',\n",
    "            'MDHm',          \n",
    "            'EX_lac__L_c']\n",
    "ID_data= [0.001674327221,\n",
    "           0.001674327221,\n",
    "           0.001674327221,\n",
    "          \n",
    "           0.0002150512027,\n",
    "           0.0002150512027,\n",
    "          0.0002150512027,\n",
    "           \n",
    "           0.002903191236,\n",
    "           0.002903191236,\n",
    "           0.002903191236,\n",
    "           0.002903191236,\n",
    "          \n",
    "           0.004178137652,\n",
    "          \n",
    "           0.003794117647,\n",
    "           0.002519171231,\n",
    "          \n",
    "           0.0006297928078,\n",
    "           0.0006297928078,\n",
    "          \n",
    "           0.002995356037,\n",
    "           0.003379376042,\n",
    "           0.003379376042,\n",
    "          0.003379376042,\n",
    "          0.003379376042,\n",
    "         0.001566801619]\n",
    "\n",
    "\n",
    "# for i in range(len(reaction)):\n",
    "#     item = core_model.reactions.get_by_id(reaction[i])\n",
    "#     a= flux_solution.loc[item]\n",
    "#     print(a)\n",
    "\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "metadata": {},
   "outputs": [],
   "source": [
    "reaction1 = ['GLCter',\n",
    "           'G6PPer',\n",
    "           'G6Pter',\n",
    "           'GLPASE1',\n",
    "           'PGMT',\n",
    "            'EX_glygn2_c',\n",
    "           'GAPD',\n",
    "           'PGK',\n",
    "           'PGM',\n",
    "            'ENO',\n",
    "            'PEPCKm',\n",
    "            'PCm',\n",
    "            'LDH_L',\n",
    "           'ME2',\n",
    "            'PYK',\n",
    "           'CSm',\n",
    "           'SUCOASm',\n",
    "            'SUCD1m',\n",
    "            'FUMm',\n",
    "            'MDHm',          \n",
    "            'EX_lac__L_c']\n",
    "\n",
    "ID_data1= [0,\n",
    "           0.00E+00,\n",
    "0,\n",
    "0,\n",
    "0,\n",
    "0,\n",
    "0,\n",
    "2,\n",
    "-2,\n",
    "-2,\n",
    "2,\n",
    "0,\n",
    "0.00E+00,\n",
    "0.00E+00,\n",
    "0.00E+00,\n",
    "2.00E+00,\n",
    "2.00E+00,\n",
    "2.00E+00,\n",
    "2.00E+00,\n",
    "2.00E+00,\n",
    "4.00E+00,\n",
    "0.00E+00]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 18,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Initial Conditions\n",
      "------------------\n",
      "GLCter: glc__D_c --> glc__D_r : (0.001674327221, 0.001674327221)\n",
      "G6PPer: g6p_r + h2o_r --> glc__D_r + pi_r : (0.001674327221, 0.001674327221)\n",
      "G6Pter: g6p_c --> g6p_r : (0.001674327221, 0.001674327221)\n",
      "GLPASE1: glygn2_c + 3.0 pi_c --> dxtrn_c + 3.0 g1p_c : (0.0002150512027, 0.0002150512027)\n",
      "PGMT: g1p_c --> g6p_c : (0.0002150512027, 0.0002150512027)\n",
      "EX_glygn2_c: glygn2_c -->  : (0.0002150512027, 0.0002150512027)\n",
      "GAPD: g3p_c + nad_c + pi_c --> 13dpg_c + h_c + nadh_c : (0.002903191236, 0.002903191236)\n",
      "PGK: 3pg_c + atp_c --> 13dpg_c + adp_c : (0.002903191236, 0.002903191236)\n",
      "PGM: 2pg_c --> 3pg_c : (0.002903191236, 0.002903191236)\n",
      "ENO: 2pg_c --> h2o_c + pep_c : (0.002903191236, 0.002903191236)\n",
      "PEPCKm: gtp_m + oaa_m --> co2_m + gdp_m + pep_m : (0.004178137652, 0.004178137652)\n",
      "PCm: atp_m + hco3_m + pyr_m --> adp_m + h_m + oaa_m + pi_m : (0.003794117647, 0.003794117647)\n",
      "LDH_L: lac__L_c + nad_c --> h_c + nadh_c + pyr_c : (0.002519171231, 0.002519171231)\n",
      "ME2: mal__L_c + nadp_c --> co2_c + nadph_c + pyr_c : (0.0006297928078, 0.0006297928078)\n",
      "PYK: adp_c + h_c + pep_c --> atp_c + pyr_c : (0.0006297928078, 0.0006297928078)\n",
      "CSm: accoa_m + h2o_m + oaa_m --> cit_m + coa_m + h_m : (0.002995356037, 0.002995356037)\n",
      "SUCOASm: atp_m + coa_m + succ_m --> adp_m + pi_m + succoa_m : (0.003379376042, 0.003379376042)\n",
      "SUCD1m: fad_m + succ_m --> fadh2_m + fum_m : (0.003379376042, 0.003379376042)\n",
      "FUMm: fum_m + h2o_m --> mal__L_m : (0.003379376042, 0.003379376042)\n",
      "MDHm: mal__L_m + nad_m --> h_m + nadh_m + oaa_m : (0.003379376042, 0.003379376042)\n",
      "EX_lac__L_c: lac__L_c -->  : (0.001566801619, 0.001566801619)\n"
     ]
    },
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>0</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0.001674</th>\n",
       "      <td>GLCter</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.001674</th>\n",
       "      <td>G6PPer</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.001674</th>\n",
       "      <td>G6Pter</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.000215</th>\n",
       "      <td>GLPASE1</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.000215</th>\n",
       "      <td>PGMT</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.000215</th>\n",
       "      <td>EX_glygn2_c</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.002903</th>\n",
       "      <td>GAPD</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.002903</th>\n",
       "      <td>PGK</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.002903</th>\n",
       "      <td>PGM</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.002903</th>\n",
       "      <td>ENO</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.004178</th>\n",
       "      <td>PEPCKm</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.003794</th>\n",
       "      <td>PCm</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.002519</th>\n",
       "      <td>LDH_L</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.000630</th>\n",
       "      <td>ME2</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.000630</th>\n",
       "      <td>PYK</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.002995</th>\n",
       "      <td>CSm</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.003379</th>\n",
       "      <td>SUCOASm</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.003379</th>\n",
       "      <td>SUCD1m</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.003379</th>\n",
       "      <td>FUMm</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.003379</th>\n",
       "      <td>MDHm</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>0.001567</th>\n",
       "      <td>EX_lac__L_c</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "                    0\n",
       "0.001674       GLCter\n",
       "0.001674       G6PPer\n",
       "0.001674       G6Pter\n",
       "0.000215      GLPASE1\n",
       "0.000215         PGMT\n",
       "0.000215  EX_glygn2_c\n",
       "0.002903         GAPD\n",
       "0.002903          PGK\n",
       "0.002903          PGM\n",
       "0.002903          ENO\n",
       "0.004178       PEPCKm\n",
       "0.003794          PCm\n",
       "0.002519        LDH_L\n",
       "0.000630          ME2\n",
       "0.000630          PYK\n",
       "0.002995          CSm\n",
       "0.003379      SUCOASm\n",
       "0.003379       SUCD1m\n",
       "0.003379         FUMm\n",
       "0.003379         MDHm\n",
       "0.001567  EX_lac__L_c"
      ]
     },
     "execution_count": 18,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "#ID_dic = [{ID_list:ID_conc}]\n",
    "print(\"Initial Conditions\\n------------------\")\n",
    "for i in range(len(reaction)):\n",
    "    item = core_model.reactions.get_by_id(reaction[i])\n",
    "    item.bounds = (ID_data[i],ID_data[i])\n",
    "    print(item, ':', item.bounds)\n",
    "    \n",
    "flux_df=pd.DataFrame(reaction,ID_data)\n",
    "flux_df"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "flux_comparison_fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))\n",
    "plot_comparison(\n",
    "    x=ID_data, y=ID_data1, compare= value,\n",
    "    observable=[reaction1], ax=ax,\n",
    "    legend=\"right outside\", plot_function=\"plot\",\n",
    "    xlim=(-16.5, 16.5), ylim=(-16.5, 16.5))\n",
    "\n",
    "flux_comparison_fig.tight_layout()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "core_model.boundary\n",
    "boundary= ['h_c',\n",
    "           'pyr_c',\n",
    "           'h2o_c',\n",
    "           'pi_c',\n",
    "           'glc__D_c',\n",
    "           'lac__L_c',\n",
    "           'co2_c',\n",
    "           'o2_c',\n",
    "           'octa_prod_c',\n",
    "           'octa_cons_c',\n",
    "           'urea_c',\n",
    "           'gln__L_c',\n",
    "           'acetone_c',\n",
    "           'bhb_c',\n",
    "           'glu__L_c',\n",
    "           'ser__L_c',\n",
    "           'cys__L_c',\n",
    "           'gly_c',\n",
    "           'glygn2_c',\n",
    "           'Tyr_ggn_c',\n",
    "           'ala__L_c',\n",
    "           'nh4_c',\n",
    "           'so3_c',\n",
    "           'glyc_3octa_c',\n",
    "           'fru_c']\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "core_model.medium\n",
    "\n",
    "medium= ['h_c', 'h2o_c','pi_c','glc__D_c','o2_c','Tyr_ggn_c','nh4_c']"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "## Looking at upper and lower bounds\n",
    "core_model.reactions.EX_glc__D_c.bounds\n",
    "\n",
    "## potentially change bounds"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "## take recon3D--> run FBA for ATPM or ATP(mitochondria)--> get fluxes--> run ensemble--> calc PERC's "
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Obtain Flux States"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "objective = \n",
    "core_model.objective = \n",
    "core_model.objective_direction = \n",
    "\n",
    "flux_solution = core_model.optimize()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#### Compare results"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# flux_comparison_fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))\n",
    "# plot_comparison(\n",
    "#     x=flux_df[\"Flux (mmol * gDW-1 * h-1)\"], y=flux_solution, compare=\"fluxes\",\n",
    "#     observable=[rid for rid in flux_df.index], ax=ax,\n",
    "#     legend=\"right outside\", plot_function=\"plot\",\n",
    "#     xlim=(-16.5, 16.5), ylim=(-16.5, 16.5),\n",
    "#     xy_line=True,\n",
    "#     xy_legend=\"best\", xlabel=\"X\", ylabel=\"Y\")\n",
    "\n",
    "# flux_comparison_fig.tight_layout()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Set paths and constants"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Define reactions"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "## want to obtain flux states, experimental data"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Indicated reactions to extract\n",
    "reaction_list = ['HEX1',\n",
    "                'PGI',\n",
    "                 'FBP',\n",
    "                'PFK',\n",
    "                'FBA',\n",
    "                'TPI',\n",
    "                'GAPD',\n",
    "                'PGK',\n",
    "                'PGM',\n",
    "                'ENO',\n",
    "                'PEPtm',\n",
    "                 'PEPCKm',\n",
    "                 'LDH_L',\n",
    "                'PYK',\n",
    "                'PCm',\n",
    "                'PYRt2m',\n",
    "                 \"ATPM\"]## Including ATPM??\n",
    "\n",
    "\n",
    "\n",
    "#                  \"DM_nadh\",\n",
    "#                  \"GLCt1\"]\n",
    "## add glucose transporter? GLCt1\n",
    "## add reactions from recon3D, \n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#### Create MASS reactions from COBRA reactions"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "new_reaction_list = []\n",
    "for item in reaction_list: \n",
    "    item = core_model.reactions.get_by_id(item)\n",
    "    new_reaction_list.append(item) \n",
    "    \n",
    "new_reaction_list\n",
    "\n",
    "# Convert cobra.Reactions to mass.MassReactions\n",
    "for rid in reaction_list:\n",
    "    reaction = core_model.reactions.get_by_id(rid)\n",
    "    glycolysis.add_reactions([MassReaction(reaction)])\n",
    "    \n",
    "## set possible exchanges for in and out?"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#### Set fluxes"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "for reaction, flux in flux_solution[glycolytic_reactions].iteritems():\n",
    "    reaction = glycolysis.reactions.get_by_id(reaction)\n",
    "#     reaction.steady_state_flux = flux * gDW_L_conversion_factor * 0.001 # mM -> M\n",
    "    print(\"{0}: {1}\".format(reaction.flux_symbol_str,\n",
    "                            reaction.steady_state_flux))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "### set bounds for both larger liver model and smaller glycolysis model\n",
    "glycolysis.reaction.EX___.lower_bound="
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Obtain Concentrations"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "ID_list = ['glc__D_c',\n",
    "           'g6p_c',\n",
    "           'f6p_c',\n",
    "           'fdp_c',\n",
    "           'dhap_c',\n",
    "           'g3p_c',\n",
    "           '13dpg_c',\n",
    "           '3pg_c',\n",
    "            '2pg_c',\n",
    "            'pep_c',\n",
    "            'pyr_c',\n",
    "            'lac__L_c',\n",
    "           'nad_c',\n",
    "            'nadh_c',\n",
    "           'amp_c',\n",
    "           'adp_c',\n",
    "            'atp_c',\n",
    "            'pi_c',\n",
    "            'h_c',          \n",
    "            'h2o_c']\n",
    "ID_conc = [10.482807,\n",
    "           0.14,\n",
    "           0.127138,\n",
    "           0.0515,\n",
    "           8.78E-03,\n",
    "           0.00878387,\n",
    "           0.176897,\n",
    "           \n",
    "           0.52063,\n",
    "           0.110561,\n",
    "           0.31,\n",
    "           0.48,\n",
    "           3.261838,\n",
    "          \n",
    "           0.0589,\n",
    "           0.0301,\n",
    "           0.365,\n",
    "           1.994952,\n",
    "            4.727146,\n",
    "           6.4,\n",
    "           0.0009,# from textbook\n",
    "           1# from textbook\n",
    "            ]\n",
    "\n",
    "#ID_dic = [{ID_list:ID_conc}]\n",
    "print(\"Initial Conditions\\n------------------\")\n",
    "for i in range(len(ID_list)):\n",
    "    item = glycolysis.metabolites.get_by_id(ID_list[i])\n",
    "    item.ic = ID_conc[i]\n",
    "    print(item, ':', item.ic)\n",
    "    \n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Load equilibrium constants\n",
    "Keq_df = pd.read_excel(\n",
    "    io=model_creation_data,\n",
    "    sheet_name=\"Equilibrium Constants\",\n",
    "    index_col=0).drop(\"Reference\", axis=1).drop(\"Stoichiometry\", axis=1)\n",
    "\n",
    "# Set equilibrium constants\n",
    "for rid, Keq in Keq_df.itertuples():\n",
    "    reaction = glycolysis.reactions.get_by_id(rid)\n",
    "    reaction.Keq = Keq\n",
    "    print(\"{0}: {1}\".format(reaction.Keq_str, Keq))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "## set upper and lower bounds\n",
    "## set constraints (constraining function)\n",
    "## set optimising function (.0*ATPM - 1.0*ATPM_reverse_5b752?? From Colton's code)\n",
    "## solve the optimal values\n",
    "## \n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#### Correct metabolite identifiers"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# for metabolite in glycolysis.metabolites:\n",
    "#     new_met_id = prefix_number_id(metabolite.id)\n",
    "#     metabolite.id = new_met_id\n",
    "# glycolysis.repair()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Set concentrations of hydrogen, water, and GDP as fixed\n",
    "for metabolite in [ \"h2o_c\", \"h_c\"]:\n",
    "    metabolite = glycolysis.metabolites.get_by_id(metabolite)\n",
    "    metabolite.fixed = True"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Keq : Set equilibrium constants"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "reaction_list = ['HEX1',\n",
    "                'PGI',\n",
    "                'PFK',\n",
    "                'FBA',\n",
    "                'TPI',\n",
    "                'GAPD',\n",
    "                'PGK',\n",
    "                'PGM',\n",
    "                'ENO',\n",
    "                 'LDH_L',\n",
    "                'PYK',\n",
    "                'ADK1',\n",
    "                 'ATPM', \n",
    "                 'DM_nadh']\n",
    "rxn_keq = [5000,\n",
    "0.32,\n",
    "1000,\n",
    "0.1,\n",
    "20.4,\n",
    "0.005334127,\n",
    "3200,\n",
    "12.7,\n",
    "4.7,\n",
    "26300,\n",
    "2220,\n",
    "1.65,\n",
    "1000000000,\n",
    "1000000000]\n",
    "\n",
    "\n",
    "\n",
    "irr_rxn=[\"SK_glc__D_c\", \"SK_amp_c\"]\n",
    "\n",
    "for i in irr_rxn:\n",
    "    item = glycolysis.reactions.get_by_id(i)\n",
    "    item.Keq=mass_config.irreversible_Keq\n",
    "    \n",
    "    \n",
    "sink=[\"SK_lac__L_c\",\"SK_pyr_c\",\"SK_h_c\",\"SK_h2o_c\"]\n",
    "# irr=[\"SK_glc__D_c\",\"SK_amp_c\"] \n",
    "for i in sink:\n",
    "    item=glycolysis.reactions.get_by_id(i)\n",
    "    item.Keq= 1\n",
    "\n",
    "\n",
    "for i in range(len(reaction_list)):\n",
    "    item = glycolysis.reactions.get_by_id(reaction_list[i])\n",
    "    item.Keq = rxn_keq[i]\n",
    "#     print(item, ':', item.Keq)\n",
    "print(\"Equilibrium Constants\\n---------------------\")\n",
    "for reaction in glycolysis.reactions:\n",
    "    print(\"{0}: {1}\".format(reaction.Keq_str, reaction.Keq))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#### Formulate QP minimization for concentrations"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "conc_solver = ConcSolver(\n",
    "    glycolysis,\n",
    "    excluded_metabolites=[\"h_c\", \"h2o_c\"],\n",
    "    constraint_buffer=1)\n",
    "\n",
    "conc_solver.setup_feasible_qp_problem(\n",
    "    fixed_conc_bounds=list(glycolysis.fixed),\n",
    "    fixed_Keq_bounds=glycolysis.reactions.list_attr(\"Keq_str\"))\n",
    "\n",
    "conc_solution = conc_solver.optimize()\n",
    "conc_solution"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#### Sample Concentrations"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# conc_solver.setup_sampling_problem(\n",
    "#     fixed_conc_bounds=list(glycolysis.fixed),\n",
    "#     fixed_Keq_bounds=glycolysis.reactions.list_attr(\"Keq_str\"))\n",
    "# for variable in conc_solver.variables:\n",
    "#     try:\n",
    "#         met = glycolysis.metabolites.get_by_id(variable.name)\n",
    "#         variable.lb, variable.ub = np.log([met.ic / 10, met.ic * 10])\n",
    "#     except:\n",
    "#         pass\n",
    "# conc_samples = sample_concentrations(conc_solver, n=n_models, seed=seed)\n",
    "# conc_samples.head()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#### Set concentrations and balance models with pseudoreactions"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# models_for_ensemble = []\n",
    "# for idx, conc_sample in conc_samples.iterrows():\n",
    "#     # Make copy of new model\n",
    "#     new_model = glycolysis.copy()\n",
    "#     new_model.id += \"_C{0:d}\".format(idx)\n",
    "#     # Get concentration sample and update model with sample\n",
    "#     new_model.update_initial_conditions(conc_sample.to_dict())\n",
    "\n",
    "#     # Determine imbalances in the reduced network\n",
    "#     fluxes = np.array(list(new_model.steady_state_fluxes.values()))\n",
    "#     imbalanced_metabolites = new_model.S.dot(fluxes)\n",
    "\n",
    "#     # Iterate through metabolites\n",
    "#     for mid, imbalance in imbalanced_metabolites.iteritems():\n",
    "#         # Ignore balanced metabolites\n",
    "#         if imbalance == 0:\n",
    "#             continue\n",
    "#         # Get metabolite object\n",
    "#         met = new_model.metabolites.get_by_id(mid)\n",
    "\n",
    "#         # Add boundary reactions for imbalanced metabolites\n",
    "#         boundary_type = \"sink\"    \n",
    "#         # Add boundary reaction with imbalance as flux value\n",
    "#         boundary_reaction = new_model.add_boundary(\n",
    "#             mid, boundary_type, boundary_condition=met.ic)\n",
    "\n",
    "#         boundary_reaction.Keq = 1\n",
    "#         if imbalance < 0:\n",
    "#             boundary_reaction.reverse_stoichiometry(inplace=True)\n",
    "#             imbalance = -imbalance\n",
    "\n",
    "#         boundary_reaction.kf = imbalance / met.ic\n",
    "#         boundary_reaction.steady_state_flux = imbalance\n",
    "#         try:\n",
    "#             # Update PERCs\n",
    "#             new_model.calculate_PERCs(\n",
    "#                 fluxes={\n",
    "#                     r: v for r, v in new_model.steady_state_fluxes.items()\n",
    "#                     if not r.boundary},\n",
    "#                 update_reactions=True)\n",
    "#         except:\n",
    "#             print(\"Negative PERCs for {0}\".format(new_model.id))\n",
    "#             continue\n",
    "#     models_for_ensemble.append(new_model)\n",
    "# print(\"Number of models in ensemble: {0:d}\".format(\n",
    "#     len(models_for_ensemble)))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Create Enzyme Modules using each concentration sample"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#### Adding from masspy documentation"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# new_reaction_list\n",
    "m=[\"adp_c\",\"amp_c\", \"atp_c\",\"pi_c\",\"nadh_c\",\"nad_c\", \"h2o_c\", \"h_c\"]\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Convert cobra.Reactions to mass.MassReactions\n",
    "for rid in m:\n",
    "    met = core_model.metabolites.get_by_id(rid)\n",
    "    glycolysis.add_metabolites([MassMetabolite(met)])"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "## Adding additional reactions\n",
    "ADK1 = MassReaction(\n",
    "    \"ADK1\",\n",
    "    name=\"Adenylate kinase\",\n",
    "    subsystem=\"Misc.\",\n",
    "    reversible=True)\n",
    "\n",
    "ADK = ['amp_c','atp_c']\n",
    "for i in ADK:\n",
    "    item = glycolysis.metabolites.get_by_id(i)\n",
    "    ADK1.add_metabolites({item:1})\n",
    "    \n",
    "ADP =['adp_c']\n",
    "for i in ADP:\n",
    "    item = glycolysis.metabolites.get_by_id(i)\n",
    "    ADK1.add_metabolites({item:-2})\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "ATPM = MassReaction(\n",
    "    \"ATPM\",\n",
    "    name=\"ATP maintenance requirement\",\n",
    "    subsystem=\"Pseudoreaction\",\n",
    "    reversible=False)\n",
    "\n",
    "ATPM_1 = ['atp_c', 'h2o_c']\n",
    "ATPM_2= ['adp_c','h_c','pi_c']\n",
    "\n",
    "for i in ATPM_1:\n",
    "    item = glycolysis.metabolites.get_by_id(i)\n",
    "    ATPM.add_metabolites({item:-1})\n",
    "    \n",
    "\n",
    "for i in ATPM_2:\n",
    "    item = glycolysis.metabolites.get_by_id(i)\n",
    "    ATPM.add_metabolites({item:1})\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "scrolled": true
   },
   "outputs": [],
   "source": [
    "DM_nadh = MassReaction(\n",
    "    \"DM_nadh\",\n",
    "    name=\"Demand NADH\",\n",
    "    subsystem=\"Pseudoreaction\",\n",
    "    reversible=False)\n",
    "\n",
    "DM_nadh_1 = ['nadh_c']\n",
    "DM_nadh_2= ['nad_c','h_c']\n",
    "\n",
    "for i in DM_nadh_1:\n",
    "    item = glycolysis.metabolites.get_by_id(i)\n",
    "    DM_nadh.add_metabolites({item:-1})\n",
    "    \n",
    "\n",
    "for i in DM_nadh_2:\n",
    "    item = glycolysis.metabolites.get_by_id(i)\n",
    "    DM_nadh.add_metabolites({item:1})\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Add new reactions\n",
    "glycolysis.add_reactions([ADK1, ATPM, DM_nadh])\n",
    "\n",
    "for reaction in glycolysis.reactions:\n",
    "    print(reaction)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Parameterize MASS model"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Boundary Reactions"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Fluxes"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "glycolysis\n",
    "glycolysis.update_S(array_type=\"DataFrame\", dtype=int)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#### Calculate PERCS"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "\n",
    "glycolysis.calculate_PERCs(update_reactions=True)\n",
    "\n",
    "print(\"Forward Rate Constants\\n----------------------\")\n",
    "for reaction in glycolysis.reactions:\n",
    "    print(\"{0}: {1:.6f}\".format(reaction.kf_str, reaction.kf))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "qcqa_model(glycolysis, parameters=True, concentrations=True,\n",
    "           fluxes=True, superfluous=True, elemental=True)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "\n",
    "# potential conversion factors?\n",
    "# enzyme modules?\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "## flux split through pathways??"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "## making enzyme modules??"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {
    "scrolled": true
   },
   "outputs": [],
   "source": [
    "## load independent fluxes\n",
    "# minspan paths, independent fluxes\n",
    "# steady state fluxes\n",
    "# print everything"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "#### Set concentrations"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Provide initial guesses for a few, setting?\n",
    "# load, set, print conc values"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "## potential optimal solutions? using QP and then comparing"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "## calculate PERC's"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "## import/merging enzyme modules into model"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Simulate at steady state"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Setup simulation object\n",
    "sim = Simulation(glycolysis, verbose=True)\n",
    "# Simulate from 0 to 1000 with 10001 points in the output\n",
    "conc_sol, flux_sol = sim.simulate(glycolysis, time=(0, 1e3, 1e4 + 1))\n",
    "# Quickly render and display time profiles\n",
    "conc_sol.view_time_profile()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Export MASS model"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "save_json_mass_model(\n",
    "    mass_model=glycolysis,\n",
    "    filename=os.path.join(model_dir, glycolysis.id + \".json\"))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "models_to_export = sim.get_model_objects(models=sim.models)\n",
    "for model in models_to_export:\n",
    "    # Save as JSON\n",
    "    save_json_mass_model(\n",
    "        mass_model=model,\n",
    "        filename=os.path.abspath(\n",
    "            os.path.join(\n",
    "                model_dir, \"JSON\", model.id + \".json\")))\n",
    "    # Save as SBML\n",
    "    write_sbml_model(\n",
    "        mass_model=model,\n",
    "        filename=os.path.abspath(\n",
    "            os.path.join(\n",
    "                model_dir, \"SBML\", model.id + \".xml.zip\")))\n",
    "# Export tables\n",
    "export_csv_files_for_models(\n",
    "    models_to_export,\n",
    "    \"for_{0}_media_{1}-{2}_isozyme_split_notebook\".format(\n",
    "        medium.lower(),\n",
    "        format_percent_str(isozyme1_percent),\n",
    "        format_percent_str(1 - isozyme1_percent)))\n",
    "print(\"Number of models exported: {0:d}\".format(len(sim.models)))"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.7.9"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
