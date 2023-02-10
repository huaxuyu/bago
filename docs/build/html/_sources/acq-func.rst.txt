Acquisition functions
---------------------

Bayesian Optimization uses acquisition function to select the next point to evaluate by real experiment. The goal of the acquisition function is to balance exploration and exploitation in the search space.

Exploration refers to the need to find a global optimum by searching different regions of the search space, 
while exploitation refers to the need to refine the search in regions where a high-performing solution is 
likely to be found. 
The acquisition function is used to decide where to sample next based on the model of the underlying function and the observed data.

There are six acquisition functions available in BAGO, including:

**Expected Improvement (EI)**: Selects the next point where the expected improvement over the current best solution is maximized.

**Probability of Improvement (PI)**: Selects the next point where the improvement over the current best solution is most probable.

**Pure exploration**: Selects the next point where the predict variance is maximized.

**Pure exploitation**: Selects the next point where the predict mean is maximized.

**Epsilon-greedy algorithm (Epsilon-greedy)**: using a parameter called "epsilon" to determine the proportion of exploratory actions.

**Upper Confidence Bound (UCB)**: Selects the next point where the balance between exploration and exploitation is optimized using a measure of uncertainty in the model.

.. note::
    
    In BAGO, Expected Improvement is used by default.

Read more about acquisition functions:

* `A Tutorial on Bayesian Optimization of Expensive Cost Functions, with Application to Active User Modeling and Hierarchical Reinforcement Learning <https://arxiv.org/abs/1012.2599>`_