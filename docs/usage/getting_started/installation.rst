Installing Kaitiaki
===================

Kaitiaki is hosted on PyPI (link here), though you can also build source from GitHub.

Prerequisites
-------------

Kaitiaki has been tested on Ubuntu 20.04 LTS. Other linux distros, and in principle any UNIX system (including Macs) should work fine. Windows is *not supported*, do not ask me to build this for a Windows machine, I will **not**.

Kaitiaki requires the Cambridge STARS code libraries. They can be downloaded `here <https://people.ast.cam.ac.uk/~stars/>`_. You will need to build the libraries into the :code:`bs` excutable. From within the STARS directory, execute the :code:`make` command.

You will need to edit :code:`run_bs`  as well. Inside the current directory, run :code:`pwd` and note the result edit the third line (first non-empty line after the shebang) to:

```
set BSDIR=<the result of your pwd command>
```

Take note of where :code:`run_bs` is. You will need it later.

The STARS code
^^^^^^^^^^^^^^

Please note that although the STARS code is *freely available*, it is not open-source as it is not released under a FOSS license. Kaitiaki is generally designed to work with the STARS code as published above. Some features however are not supported. Kaitiaki is capable of modifying STARS in-situ and can do this on request to enable such features. Perhaps one day we will be able to release the libraries alongside Kaitiaki in a FOSS-compatible license.

From Source
-----------

Building from source is straightforward:

1. Fork the `Github Repo <https://github.com/Krytic/kaitiaki>`_
2. Locally, run :code:`git clone <your url>`
3. In a terminal, navigate to the directory you installed kaitiaki to
4. Run :code:`pip install -e .` to install it as an editable PIP file

From PyPI
---------

Kaitiaki is currently *not* on PyPI (yet), but when it is released you will be able to use :code:`pip install kaitiaki` to install it.