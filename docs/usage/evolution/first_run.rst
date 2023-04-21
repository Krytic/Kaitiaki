Your First Kaitiaki Run
=======================

The simplest example of Kaitiaki is to run the evolution of a :math:`1 M_\odot`, solar metallicity star. This can be done with the following snippet:

.. code-block::

	import kaitiaki

	STARS = kaitiaki.STARS.STARSController()
	STARS.blit()
	STARS.load_default_modin()
	STARS.run()

Let's examine these in more detail.

How STARS works
---------------

STARS requires fundamentally three customisable things:

1. A model input
2. A data file
3. Opacity tables.

modin
^^^^^

The model input -- often shortened to the :code:`modin` file -- is the state of the model upon start of simulation. Kaitiaki ships a sample model input -- for :math:`1~M_\odot,~Z=0.020` (if you want to do something more interesting than the sun, you need to inflate it yourself -- see a future tutorial). The :code:`STARS.load_default_modin()` creates that file on disk for you.

data
^^^^

The data file is a list of potential inputs. The STARS manual has a more detailed list of parameters, but the default ones are set up to allow you to just immediately evolve a :math:`1~M_\odot` model until roughly the AGB phase. This is loaded in :code:`STARS.blit()`.

COtables
^^^^^^^^

The opacity tables are what enable you to run different metallicities. For instance, if you wanted to run a :math:`Z=0.014` model, you would need to load the z014 COtable. This is loaded in :code:`STARS.blit()`.