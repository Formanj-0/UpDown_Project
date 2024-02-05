import scr.UpDown_Module.preprocessing_functions as ppf
import os
import scanpy as sc
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import ast 




class overlord_pp():
    def __init__(self) -> None:
        self.adata = None
        self.data_path = None
        self.data_type = None
        self.steps = []
        self.params = {}



    # Decorator to record method calls and their parameters
    def record_method_calls(func):
        def wrapper(*args, **kwargs):
            # Record the method name
            method_name = func.__name__
            overlord = args[0]  # Assuming the first argument is the overlord_pp instance
            overlord.steps.append((method_name))  # Append method name to self.steps

            # Record the passed variables
            overlord.params[method_name] = kwargs  # Update self.params with method parameters

            # Call the method
            result = func(*args, **kwargs)

            return result

        return wrapper


    @record_method_calls
    def load_data(self, data_path, data_type):
        """
        Loads data from the specified data_path and data_type using the preprocessing_functions module.
        Updates the self.data_path, self.data_type, and self.adata attributes.
        Returns the loaded data.
        """
        self.data_path = data_path
        self.data_type = data_type
        self.adata = ppf.load_data(data_path, data_type)
        return self.adata


    @record_method_calls
    def save_data(self, data_path, data_type, foldername, cwd):
        """
        Saves the self.adata object to the specified data_path.
        Updates the self.data_path and self.data_type attributes.
        Writes self.steps and self.params to separate files.
        Returns the saved data.
        """
        self.data_path = data_path
        self.data_type = data_type
        self.adata.write(data_path)

        # Create the settings folder if it does not exist
        settings_folder = os.path.join(cwd, f"settings_{foldername}")
        if not os.path.exists(settings_folder):
            os.makedirs(settings_folder)

        # Save self.steps to a file
        steps_file_path = os.path.join(settings_folder, "steps.txt")
        with open(steps_file_path, "w") as steps_file:
            steps_file.write(str(self.steps))

        # Save self.params to a file
        params_file_path = os.path.join(settings_folder, "params.txt")
        with open(params_file_path, "w") as params_file:
            params_file.write(str(self.params))

        return self.adata


    @record_method_calls
    def load_and_run_files(self, cwd, foldername):
        """
        Loads and runs the files located in the specified location.
        Uses the functions in the `steps` list and the parameters from the `params` dictionary.
        """
        # Define the file paths
        steps_file_path = os.path.join(cwd, f"settings_{foldername}", "steps.txt")
        params_file_path = os.path.join(cwd, f"settings_{foldername}", "params.txt")

        # Load the steps
        with open(steps_file_path, "r") as steps_file:
            steps_str = steps_file.read()
        steps = ast.literal_eval(steps_str)

        # Load the params
        with open(params_file_path, "r") as params_file:
            params_str = params_file.read()
        params = ast.literal_eval(params_str)

        print(f"Loaded steps: {steps}")
        print(f"Loaded params: {params}")
        # Run the functions in order
        for step in steps:
            method_name = step
            method_params = params.get(method_name, {})
            print(f"Running {method_name} with parameters {method_params}")
            getattr(self, method_name)(**method_params)



    @record_method_calls
    def remove_cells(self, min_genes=None, max_genes=None, min_counts=None, max_counts=None):
        """
        Removes cells from the self.adata object based on the specified gene and count thresholds.
        Updates the self.adata object.
        """
        print("===== Removing cells =====")
        cells_intial = self.adata.n_obs
        genes_initial = self.adata.n_vars
        self.adata = ppf.remove_cells(self.adata, min_genes, max_genes, min_counts, max_counts)
        cells_final = self.adata.n_obs
        genes_final = self.adata.n_vars
        print(f"Removed {cells_intial - cells_final} cells")
        print(f"Removed {genes_initial - genes_final} genes")
        return self.adata



    @record_method_calls
    def remove_genes(self, min_cells=None, max_cells=None, min_counts=None, max_counts=None):
        """
        Removes genes from the self.adata object based on the specified cell and count thresholds.
        Updates the self.adata object.
        """
        print("===== Removing genes =====")
        genes_initial = self.adata.n_vars
        cells_initial = self.adata.n_obs
        self.adata = ppf.remove_genes(self.adata, min_cells, max_cells, min_counts, max_counts)
        genes_final = self.adata.n_vars
        cells_final = self.adata.n_obs
        print(f"Removed {genes_initial - genes_final} genes")
        print(f"Removed {cells_initial - cells_final} cells")
        return self.adata



    @record_method_calls
    def label_qc_genes(self, label_list, search_list):
        """
        Labels genes in the self.adata.var dataframe based on the specified label_list and search_list.
        Updates the self.adata.var dataframe.
        """
        for label, search_term in zip(label_list, search_list):
            self.adata.var[label] = self.adata.var_names.str.startswith(search_term)



    @record_method_calls
    def calculate_qc_metrics(self, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True):
        """
        Calculates QC metrics for the self.adata object.
        """
        old_qc_obs_names = list(self.adata.obs.columns)  # Get the current QC observation names
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=qc_vars, percent_top=percent_top, log1p=log1p, inplace=inplace)
        new_qc_obs_names = list(set(self.adata.obs.columns) - set(old_qc_obs_names))  # Get the new QC observation names
        self.adata.uns['qc_obs_names'] = new_qc_obs_names  # Store the new QC observation names in adata.uns['qc_obs_names']



    @record_method_calls
    def visualize_qc_metrics(adata, method='standard'):
        if method == 'standard':
            sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
            sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
            sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
        elif method == 'Jacks':
            import matplotlib.pyplot as plt
            plt.scatter(adata.obs['n_genes_by_counts'], adata.obs['total_counts'], c=adata.obs['pct_counts_mt'])
            plt.xlabel('genes by counts')
            plt.ylabel('total number of reads')
            plt.title(adata.uns['Dataset_name'] + ' before QC')
            plt.colorbar()
            plt.show()
        else:
            print("Invalid visualization method. Please choose a valid method.")


    @record_method_calls
    def simple_qc(self, method, value, column_name, operator):
        if method == 'comparison':
            self.adata = ppf.remove_cells_by_qc_var(self.adata, column_name, value, operator)
        # Add future methods here



    @record_method_calls
    def normalize_counts(self, target_sum=1e4):
        """
        Normalizes the counts using sc.pp.normalize_total.
        Allows the target sum to be changed.
        Plots histograms of the sums before and after normalization in both normal linear space and log space.
        """
        # Calculate sums before normalization
        sums_before = self.adata.X.sum(axis=1)

        # Normalize counts
        sc.pp.normalize_total(self.adata, target_sum=target_sum)

        # Calculate sums after normalization
        sums_after = self.adata.X.sum(axis=1)

        # Plot histograms
        fig, axs = plt.subplots(2, 2, figsize=(10, 8))

        # Plot histogram of sums before normalization in normal linear space
        axs[0, 0].hist(sums_before, bins=50)
        axs[0, 0].set_title('Sums Before Normalization (Linear Space)')
        axs[0, 0].set_xlabel('Sum')
        axs[0, 0].set_ylabel('Frequency')

        # Plot histogram of sums before normalization in log space
        axs[0, 1].hist(sums_before, bins=50)
        axs[0, 1].set_title('Sums Before Normalization (Log Space)')
        axs[0, 1].set_xlabel('Log(Sum)')
        axs[0, 1].set_ylabel('Frequency')
        axs[0, 1].set_yscale('log')
        axs[0, 1].set_xscale('log')

        # Plot histogram of sums after normalization in normal linear space
        axs[1, 0].hist(sums_after, bins=50)
        axs[1, 0].set_title('Sums After Normalization (Linear Space)')
        axs[1, 0].set_xlabel('Sum')
        axs[1, 0].set_ylabel('Frequency')

        # Plot histogram of sums after normalization in log space
        axs[1, 1].hist(sums_after, bins=50)
        axs[1, 1].set_title('Sums After Normalization (Log Space)')
        axs[1, 1].set_xlabel('Log(Sum)')
        axs[1, 1].set_ylabel('Frequency')
        axs[1, 1].set_yscale('log')
        axs[1, 1].set_xscale('log')

        plt.tight_layout()
        plt.show()



    @record_method_calls
    def convert_to_log1p(self):
        """
        Saves the current counts to adata.layers["Counts"] and converts to log1p space using sc.pp.log1p.
        """
        # Save current counts to adata.layers["Counts"]
        self.adata.layers["Counts"] = self.adata.X.copy()

        # Convert to log1p space
        sc.pp.log1p(self.adata)


    @record_method_calls
    def move_column_to_obs(self, column_name):
        """
        Moves a column from adata.X array to a column with the same name in adata.obs dataframe,
        and removes the column from adata.X array.
        
        Args:
            adata (AnnData): Annotated data object.
            column_name (str): Name of the column to be moved.
        """
        # Check if column_name exists in adata.X
        if column_name in self.adata.var_names:
            # Get the index of the column in adata.X
            column_index = self.adata.var_names.get_loc(column_name)
            
            # Get the column data from adata.X
            column_data = self.adata.X[:, column_index]
            
            # Add the column to adata.obs with the same name
            self.adata.obs[column_name] = column_data
            
            # Remove the column from adata.X
            self.adata = self.adata[:, ~column_index]
        else:
            print(f"Column '{column_name}' does not exist in adata.X")


    @record_method_calls
    def move_row_to_var(self, row_name):
        """
        Moves a row from adata.obs dataframe to a row with the same name in adata.var dataframe,
        and removes the row from adata.obs dataframe.
        
        Args:
            adata (AnnData): Annotated data object.
            row_name (str): Name of the row to be moved.
        """
        # Check if row_name exists in adata.obs
        if row_name in self.adata.obs_names:
            # Get the index of the row in adata.obs
            row_index = self.adata.obs_names.get_loc(row_name)
            
            # Get the row data from adata.obs
            row_data = self.adata.obs.loc[row_name]
            
            # Add the row to adata.var with the same name
            self.adata.var[row_name] = row_data
            
            # Remove the row from adata.obs
            self.adata = self.adata[~row_index, :]
        else:
            print(f"Row '{row_name}' does not exist in adata.obs")























