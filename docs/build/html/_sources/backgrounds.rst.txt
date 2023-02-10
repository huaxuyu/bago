Backgrounds
-----------

Bayesian optimization is a probabilistic model-based optimization algorithm 
for finding the maximum or minimum of an unknown objective function. It's an 
iterative method that uses Bayesian inference to model the underlying function 
based on the observations made during the optimization process.

The algorithm starts with an initial set of points and builds a statistical 
model of the objective function based on the observed data. The algorithm 
then decides the next point to evaluate based on the model and an acquisition 
function, which balances exploration (sampling in regions where the model is uncertain) 
and exploitation (sampling in regions where the model predicts high values). 
This process is repeated until the algorithm reaches a stopping criterion or a maximum number of iterations.

Bayesian optimization is particularly well-suited to optimization problems with 
expensive objective functions, such as simulations, experiments, or real-world 
applications, where evaluating the function is time-consuming or costly. 
It has been used in various fields, including machine learning, computer vision, 
robotics, and engineering, to optimize hyperparameters of models, design experiments, 
and find optimal control policies. In our case, evaluating an elution gradient 
using LC-MS is a time-consuming process, which makes Bayesian optimization a 
great tool for optimizing the gradient.

Read more about Bayesian optimization:

* `A Tutorial on Bayesian Optimization of Expensive Cost Functions, with Application to Active User Modeling and Hierarchical Reinforcement Learning <https://arxiv.org/abs/1012.2599>`_