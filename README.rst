
AuFAMe : Automated Functional Annotation and Metabolic networks building from bacterial genomes. 
==========================================

Workflow to reconstruct multiple metabolic graphs directly from sequence fasta.

.. contents:: Table of contents
   :backlinks: top
   :local:

License
--------
This workflow is licensed under the GNU GPL-3.0-or-later, see the `LICENSE <https://github.com/AuReMe/aufame/blob/main/LICENSE>`__ file for details.

Installation
------------

Dependencies
~~~~~~~~~~~~

These tools are needed:

	- `Prokka <https://github.com/tseemann/prokka>`__

	- `EggNOG-mapper <https://github.com/eggnogdb/eggnog-mapper>`__

	- `Bakta <https://github.com/oschwengers/bakta>`__

	- `Pathway Tools <http://bioinformatics.ai.sri.com/ptools/>`__ (which needs `Blast <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`__)


And some python packages:

	- `mpwt <https://github.com/AuReMe/mpwt>`__

	- `padmet <https://github.com/AuReMe/padmet>`__

	- `pandas <https://pandas.pydata.org/>`__


Installation of Pathway Tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run annotation based reconstruction, you need to install Pathway Tools. This tool is 
available at the `Pathway Tools <http://bioinformatics.ai.sri.com/ptools/>`__ website. 


Getting the MetaCyc PADMET file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You also should install the MetaCyc_XX.X.padmet (the version number of 
`MetaCyc <https://metacyc.org/>`__  is replaced with XX.X), and then you should update your 
config.txt files for each study. This is the way to get a MetaCyc_XX.padmet file: 
Firstly, download the flat files of `MetaCyc <https://metacyc.org/>`__ in DAT format at the
`https://biocyc.org/download.shtml <https://biocyc.org/download.shtml>`__ webpage. 
Secondly, put all the downloaded DAT files in a directory (it is named FLAT_DIR here). 
Thirdly run this command:

.. code:: sh

	padmet pgdb_to_padmet --pgdb=FLAT_DIR --output=metacyc_XX.X.padmet --version=XX.X --db=Metacyc -v

pip
~~~

If you have installed all the dependencies, you can just install AuFAMe with:

.. code:: sh

	pip install aufame

Usage
-----

Before launching AuFAMe, you will need to setup the configfile and place it in the output directory of your choice. 
Configfile is composed of 5 categories : 

	- **Conda environments :** provides path to yaml files (by default) or already available conda environment.
	
	- **Databases :** provides path to the 3 databases used during the workflow : Bakta, EggNOG-mapper, and MetaCyc.
	
	- **Files :** provides path to :
		- input directory containing genomes ; each genome must be placed in a subdirectory having the same name 
		- taxon file ; it can follow `Prolipipe's taxonfile structure <https://github.com/AuReMe/prolipipe?tab=readme-ov-file#preparation-of-taxfile>`__, the 3 columns are mandatory (but the first, "Species", can be empty)
		- output directory : not to change if the config file is already in your output directory
	
	- **EggNOG optimization option :** : AuFAMe can place temporarily the EggNOG database in memory to accelerate annotation. It's composed of 2 parameters : 
		- `eggnog_db_optim` : fill it with `"yes"` if you want AuFAMe to move it temporarily in memory (needs 50G free).
		- `eggnog_db_optim_path` : you can specify the directory name in memory (need to start with /dev/shm)
	
	- **Other parameters :** : you can specify the annotation tools to rely on, and intermediate files to remove 

Once in an environment able to run Snakemake, you can launch AuFAMe with the following command :

.. code:: bash

	snakemake --use-conda --conda-frontend conda -s path/to/aufame.smk -d [OUTPUT_DIR] -p --rerun-triggers mtime --resources process_data_jobs=10 --rerun-incomplete --keep-going --cores [NB_THREADS]