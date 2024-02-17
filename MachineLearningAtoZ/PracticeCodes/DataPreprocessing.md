Importing and Prepocessing

```python
# Libraries required to be imported
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split

# Iris dataset will be example dataset to play here. Load it as a dataframe using pandas (pd)
# More info: https://scikit-learn.org/stable/auto_examples/datasets/plot_iris_dataset.html
dataset = pd.read_csv("iris.csv")

# Let's see the first five rows of the dataset
dataset.head()

# Features will be created in the matrix format and variable in the vector format.
# iloc allows to subset the relevant columns and matrix
# [:,:-1] all rows and all columns except the last column
# [:, -1] al rows of the last column only
X = dataset.iloc[: , :-1].values #features matrix
y = dataset.iloc[: , -1].values #dependent variable vector

# Print the X and y
print(X)
print(y)
```

Handle Missing Data
