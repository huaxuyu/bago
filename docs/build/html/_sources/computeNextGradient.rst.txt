Compute next gradient
---------------------

A function in the MSData object to find the top signals in the MS1. 

.. function:: computeNextGradient(self, acqFunc="ei")

    Function to calculate the next gradient to run using an acquisition function.

    :param acqFunc: str. Acquisition functions in lower case. Valid options are: ei (expected improvement), explore (pure exploration), exploit (pure exploitation), eps (epsilon-greedy algorithm), pi (probability of improvement), rand (random a gradient), ucb (upper confidence bound)

To use this function:

.. code-block::

    # You need a gpModel object (model)

    model.computeNextGradient()


