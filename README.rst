MASS Liver
==========
This repository is for the MASSpy liver model project.

Clone the repository using git clone, then navigate to the cloned repository and run the following command::

    pip install -e ".[all]"

Afterwards, make sure all worked by running the following from the repository::

    tox

Development Notes
-----------------
* Sicily: Glycolysis, Gluconeogenesis, PPP
* Shuling: Sucrose Metabolism, Fructolysis, Glycogensis, Glycogenolysis
* Riya: Pyruvate Dehydrogenase, TCA Cycle, Oxidative Phosphorylation, Ketogenesis, Ketone Oxidation
* Mahima: Amino Acid Degradation, Amino Acid Synthesis, Urea Cycle, Aspartate- Malate Shuttle
* TBD: Citrate Malate Shuttle


Dev. Branches
+++++++++++++
* main: Managed by Zack, main branch for stable code
* devel: Managed by Zack, main branch for development of the liver model

* Create a branch called ``devel-*`` where the ``*`` represents whatever you want to call the branch. Make sure this branch is from the ``devel`` branch. 

    1. git checkout devel
    2. git branch devel-zh-patch
    3. git checkout devel-zh-patch

* When you are ready to merge the branch into ``devel``, create a pull request and request a review from ``z-haiman``.

To speed up the process of getting your code approved and working...

* Before pushing local changes to the GitHub repository, try testing them locally using tox::

    tox

* Sometimes you may need to reset the tox cache. To do this, run ``tox --recreate`` and all the environment caches will be rebuilt.
