==============
python-scripts
==============

Python scripts for hydrodynamic simulations

Process description
==================
Create Python virtual environment::

python3 -m venv --clear venv
./venv/bin/pip install --upgrade pip
./venv/bin/pip install -r requirements.txt

To run scripts, use for example::

./venv/bin/spectrum.py "directory" "output" NoF NoE

Contribution guidelines
=======================

This repository uses Ruff_ to ensure consistent code style and for linting.
Shell scripts use shfmt_ and other supported file types use Prettier_ for
formatting. Make sure to run ::

   ruff check --fix .
   ruff format .
   shfmt -l -w .
   prettier --write --ignore-unknown "**/*"

Alternatively, set up pre-commit_ to take care of these tasks. First, install
pre-commit::

   pip install pre-commit

Then, after cloning the repository, install the Git hook scripts::

   pre-commit install

.. _Ruff: https://github.com/astral-sh/ruff
.. _shfmt: https://github.com/mvdan/sh
.. _Prettier: https://github.com/prettier/prettier
.. _pre-commit: https://pre-commit.com
