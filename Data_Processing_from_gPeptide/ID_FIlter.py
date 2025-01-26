
import pandas as pd


def filter_expression_matrix(file_name: str, coverage_threshold: int) -> None:
    data = pd.read_csv(file_name, index_col=0)
    samples = list(data.columns)
    selected_gpids_data: dict[str, pd.Series] = {}
    for gpid in data.index:
        if gpid in selected_gpids_data.keys():
            raise Exception(f'Repetitive ID detected: {gpid}')
        line_data = data.loc[gpid]
        count = 0
        for sample in samples:
            if line_data[sample] > 0.0: count += 1
        if count > coverage_threshold:
            selected_gpids_data[gpid] = line_data
    selected_data: pd.DataFrame = pd.DataFrame(columns=samples)
    for gpid in selected_gpids_data.keys():
        selected_data.loc[gpid] = selected_gpids_data[gpid].values
    selected_data.to_csv('Expression_Matrix.csv', index=True)
    
    
if __name__ == '__main__':
    
    filter_expression_matrix('test_expression_matrix.csv', 10)