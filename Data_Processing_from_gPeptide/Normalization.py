
import numpy as np
import pandas as pd


def minmax_scale(data: pd.DataFrame):
    total_values = []
    for col in data.columns: total_values.extend(data[col].values)
    min_value = min(total_values)
    max_value = max(total_values)
    base = max_value - min_value
    for c in data.columns:
        col_data = data[c]
        data[c] = [(x - min_value) / base for x in col_data.values]
        

def sigmoid_scale(data: pd.DataFrame, cons: float = 3.5) -> None:
    total_values = []
    for col in data.columns: total_values.extend(data[col].values)
    median = np.percentile(total_values, 50)
    for c in data.columns:
        data[c] = 1 / (1 + np.exp(- cons / abs(median) * (data[c].values - median)))


if __name__ == '__main__':
    
    expression_matrix_imputed = pd.read_csv('Expression_Matrix_BEC.csv', index_col=0)
    minmax_scale(expression_matrix_imputed)
    expression_matrix_imputed.to_csv('Expression_Matrix_Norm.csv')
    sigmoid_scale(expression_matrix_imputed, cons=1.0)
    expression_matrix_imputed.to_csv('Expression_Matrix_for_Cluster.csv')