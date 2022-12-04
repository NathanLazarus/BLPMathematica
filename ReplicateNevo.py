import pyblp
import numpy as np
import pandas as pd

pyblp.options.digits = 2
pyblp.options.verbose = False
print(pyblp.__version__)

product_data = pd.read_csv(pyblp.data.NEVO_PRODUCTS_LOCATION)
print(product_data.head())

agent_data = pd.read_csv(pyblp.data.NEVO_AGENTS_LOCATION)
print(agent_data.head())

agent_formulation = pyblp.Formulation('0 + income + income_squared + age + child')
print(agent_formulation)

X1_formulation = pyblp.Formulation('0 + prices', absorb='C(product_ids)')
X2_formulation = pyblp.Formulation('1 + prices + sugar + mushy')
product_formulations = (X1_formulation, X2_formulation)
print(product_formulations)

# mc_integration = pyblp.Integration('monte_carlo', size=50, specification_options={'seed': 0})
# print(mc_integration)

# pr_integration = pyblp.Integration('product', size=5)
# print(pr_integration)


initial_sigma = np.diag([0.3302, 2.4526, 0.0163, 0.2441])
initial_pi = np.array([
  [ 5.4819,  0,      0.2037,  0     ],
  [15.8935, -1.2000, 0,       2.6342],
  [-0.2506,  0,      0.0511,  0     ],
  [ 1.2650,  0,     -0.8091,  0     ]
])

# print('product_data.shape', product_data.shape)
# print('agent_data.shape', agent_data.shape)

X1_formulation_with_d = pyblp.Formulation('0 + prices + C(product_ids)')
product_formulations_with_d = (X1_formulation_with_d, X2_formulation)
nevo_problem_with_d = pyblp.Problem(
    product_formulations_with_d,
    product_data,
    agent_formulation,
    agent_data
)
nevo_results_with_d = nevo_problem_with_d.solve(
    initial_sigma,
    initial_pi,
    optimization=pyblp.Optimization('bfgs', {'gtol': 0.5}),  # use a loose tolerance as in the original paper
    method='1s'
)

# print(nevo_results_with_d)
# # collect the estimated fixed effects, their covariances, and the desired regressors
# print('nevo_results_with_d.__dict__', nevo_results_with_d.__dict__)
# print('nevo_results_with_d.parameter_covariances', nevo_results_with_d.parameter_covariances)
# print('nevo_results_with_d.parameter_covariances.shape', nevo_results_with_d.parameter_covariances.shape)
# print('nevo_results_with_d.beta', nevo_results_with_d.beta)
d = nevo_results_with_d.beta[1:]
V = nevo_results_with_d.parameter_covariances[-d.size:, -d.size:]
X = np.zeros((d.size, 3))
for i, label in enumerate(nevo_results_with_d.beta_labels[1:]):
    row = product_data[product_data['product_ids'] == label.split("'")[1]].iloc[0]
    X[i] = np.r_[1, row[['sugar', 'mushy']]]

# print('X', X)
# form the minimum distance estimates and their standard errors
W = np.linalg.inv(V)
beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ d)
beta_se = np.sqrt(np.diag(np.linalg.inv(X.T @ W @ X)) / nevo_problem_with_d.N)
print(np.c_[beta, beta_se])