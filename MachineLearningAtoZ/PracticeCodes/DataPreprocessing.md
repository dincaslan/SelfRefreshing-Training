**Importing and Prepocessing**

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

**Handling Missing Data**

```python
# Libraries required to be imported
from sklearn.impute import SimpleImputer
import pandas as pd
import numpy as np

# Load the dataset as dataframe using pandas as pd
# Similar example from another tutorial: https://monashdatafluency.github.io/python-workshop-base/modules/missing_values/
# In our case, we have a dataset having 10 columns and more than 700 rows. The first 9 columns are numerical values consisting the risk factors of the diabetes patient,the last columns is the binary data with the outcome as 0 or 1 (either present or not)
newdata = pd.read_csv("diabetes.csv")
X = newdata.iloc[: , :-1].values #features matrix
y = newdata.iloc[: , -1].values #outcome

# Identify the missing data which is represented as NaN
# pandas isnull() function helps to recognize the missing values
newdata[pd.isnull(newdata).any(axis=1)]

# Printing the number of missing entries in each column
print(len(newdata[pd.isnull(newdata)]))

# Configure the SimpleImputer class, fit the imputer on the dataframe, and finally transform the dataframe with updated version of the missing values
# You can drop the NAs, dataset.dropna()
# Mean, median or constant value options might be choosen to replace the missing value in the simple imputer instance
# This step requires numpy to be imported as np
imputer = SimpleImputer(missing_values=np.nan, strategy="mean")
imputer.fit(X[:, 1:9])
X[:, 1:9] = imputer.transform(X[:, 1:9])

#Print your updated features
print(X)
```

 **Encoding Categorical Data**

```python
# We will use OneHotEncoder to make categorical data numerical. In this case, it would be better to be aware of advantages and disadvantages of the method. If you want to learn more: https://www.geeksforgeeks.org/ml-one-hot-encoding-of-datasets-in-python/

# Libraries required to be imported
import pandas as pd
import numpy as np

```
 
