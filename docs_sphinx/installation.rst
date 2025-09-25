Installation Guide
==================

This guide will walk you through the installation process for ClearWater Modules.

Requirements
------------

ClearWater-modules was developed with **Python 3.11** and requires:

- Python 3.11 or higher
- Conda or Miniconda for environment management
- Git (for cloning the repository)

Step-by-Step Installation
-------------------------

1. Install Miniconda or Anaconda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We recommend installing the lightweight `Miniconda <https://docs.conda.io/projects/miniconda/en/latest/>`_ that includes Python, the conda environment and package management system, and their dependencies.

.. note::
   Follow conda defaults to install in your local user directory. **DO NOT** install for all users, to avoid substantial headaches with permissions.

If you have already installed the `Anaconda Distribution <https://www.anaconda.com/download>`_, you can use it to complete the next steps, but you may need to update to the latest version.

2. Clone or Download the Repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Clone the repository using Git:

.. code-block:: bash

   git clone https://github.com/EcohydrologyTeam/ClearWater-modules.git
   cd ClearWater-modules

Alternatively, download the ZIP file from the GitHub repository and extract it to your desired location.

3. Create a Conda Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a dedicated environment using the provided ``environment.yml`` file:

.. code-block:: bash

   conda env create --file environment.yml

For faster installation, use the libmamba solver:

.. code-block:: bash

   conda env create -f environment.yml --solver=libmamba

Activate the environment:

.. code-block:: bash

   conda activate clearwater-modules

4. Add ClearWater-modules to Python Path
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To make the ``clearwater_modules`` package accessible in your Python environment, use conda develop:

.. code-block:: bash

   conda develop /path/to/ClearWater-modules

Replace ``/path/to/ClearWater-modules`` with the actual path to your cloned repository.

Updating the Environment
------------------------

To update your environment with the latest dependencies:

.. code-block:: bash

   conda env update -f environment.yml --solver=libmamba --prune

To recreate the environment from scratch:

.. code-block:: bash

   conda env create -f environment.yml --solver=libmamba --force

Verifying Installation
----------------------

To verify that the installation was successful, activate your environment and try importing the modules:

.. code-block:: python

   import clearwater_modules
   from clearwater_modules import tsm, nsm1, nsm2, gsm, csm, msm

If no errors occur, the installation is complete!

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**Import Errors**: Ensure that you've run ``conda develop`` with the correct path to the repository.

**Environment Activation**: Make sure to activate the conda environment before using the modules:

.. code-block:: bash

   conda activate clearwater-modules

**Permission Errors**: If you encounter permission errors, ensure you installed conda in your user directory, not system-wide.

Getting Help
~~~~~~~~~~~~

If you encounter issues not covered here:

1. Check the GitHub issues page
2. Create a new issue with details about your problem
3. Include your system information and error messages