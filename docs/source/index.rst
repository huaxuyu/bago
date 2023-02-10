.. BAGO documentation master file, created by
   sphinx-quickstart on Tue Feb  7 13:47:18 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to BAGO's documentation!
================================

`BAGO`_ is a Python package for Bayesian Optimization of Liquid Chromatographic Elution Gradient.
Use BAGO to design a gradient for your LC-MS/MS analysis today!

BAGO enables:

Highly efficient gradient optimization
    Find an optimal gradient for your LC-MS/MS analysis within 10 runs.
    Wonder why BAGO is efficient? Read more about :doc:`/acq-func`.

Omics-scale evaluation on compound separation
    Separation efficiency was defined to evaluate the performance of a gradient.
    Wonder how omics-scale evaluation is achieved? Read more about :doc:`/encodings`.

Broader discovery of chemical space
    Expand your discovery of chemical space by improving identification and quantification.
    Wonder how BAGO can help you? Read more about :doc:`/applications`.

.. _BAGO: https://github.com/Waddlessss/BAGO


Get Started
-----------

Start your journey with BAGO by reading the following pages:

* **Background**: 
  :doc:`/backgrounds`

* **Getting Started**:
  :doc:`With A Jupyter Notebook </getting-started-with-ipynb>` |
  :doc:`With A GUI Software </getting-started-with-software>` |

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Get Started

   /backgrounds
   /getting-started-with-ipynb
   /getting-started-with-software
   /encodings
   /applications


BAGO Functions
--------------

Learn more about the functions in BAGO.

* **Functions to manipulate MS data**:
  :doc:`/MSdata-obj` |
  :doc:`/readMSdata` |
  :doc:`/computeSecondGradient` |
  :doc:`/findTopSignals` |
  :doc:`/computeSepEff` |
  :doc:`/getBPCData` |
  :doc:`/plotBPC` |
  :doc:`/dotProd` |
  :doc:`/numberUniqueMS2` |
  :doc:`/getUniqueMz` |
  :doc:`/getMobilePhasePct` |
  :doc:`/outputConfig`

* **Functions to build Bayesian optimization model**:
  :doc:`/gpModel-obj` |
  :doc:`/fit` |
  :doc:`/genSearchSpace` |
  :doc:`/computeNextGradient` |
  :doc:`/updateModel` |
  :doc:`/acq-func`

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: BAGO Functions

   /MSdata-obj
   /readMSdata
   /computeSecondGradient
   /findTopSignals
   /computeSepEff
   /getBPCData
   /plotBPC
   /dotProd
   /numberUniqueMS2
   /getUniqueMz
   /getMobilePhasePct
   /outputConfig
   /gpModel-obj
   /fit
   /genSearchSpace
   /computeNextGradient
   /updateModel
   /acq-func


Useful Links
------------

* **BAGO on GitHub**:
   `BAGO source code and more details <https://github.com/Waddlessss/BAGO>`_

* **BAGO on PyPI**:
   `BAGO Python package <https://pypi.org/project/bago/>`_


Citation
--------

Please cite:

* **BAGO paper**:


Further Reading
---------------

* **Bayesian optimization**:

   - `A Tutorial on Bayesian Optimization <https://arxiv.org/abs/1807.02811>`_
   - `A Tutorial on Bayesian Optimization of Expensive Cost Functions, with Application to Active User Modeling and Hierarchical Reinforcement Learning <https://arxiv.org/abs/1012.2599>`_

* **Gaussian process regression**:

   - `Gaussian process in scikit-learn <https://scikit-learn.org/stable/modules/gaussian_process.html>`_
   - `Gaussian Processes for Machine Learning <https://gaussianprocess.org/gpml/>`_

* **Liquid chromtography**:

   - `Introduction to Modern Liquid Chromatography <https://onlinelibrary.wiley.com/doi/book/10.1002/9780470508183>`_