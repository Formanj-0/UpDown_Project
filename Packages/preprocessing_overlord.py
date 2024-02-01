import preprocessing_functions as ppf
import os
import scanpy as sc
import numpy as np
import anndata as ad
import ast 
import os


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

            # Record the passed variables
            overlord.params[method_name] = kwargs

            # Call the method
            result = func(*args, **kwargs)

            return result

        return wrapper



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



    def save_data(self, data_path, data_type):
        """
        Saves the self.adata object to the specified data_path.
        Updates the self.data_path and self.data_type attributes.
        Writes self.steps and self.params to separate files.
        Returns the saved data.
        """
        self.data_path = data_path
        self.data_type = data_type
        self.adata.write(data_path)

        # Save self.steps to a file
        steps_file_path = f"settings_{os.path.splitext(__file__)[0]}"
        with open(steps_file_path, "w") as steps_file:
            steps_file.write(str(self.steps))

        # Save self.params to a file
        params_file_path = f"settings_{os.path.splitext(__file__)[0]}"
        with open(params_file_path, "w") as params_file:
            params_file.write(str(self.params))

        return self.adata



    def load_and_run_files(self):
        """
        Loads and runs the files located in the specified location.
        Uses the functions in the `steps` list and the parameters from the `params` dictionary.
        """
        # Load the steps from the file
        steps_file_path = f"settings_{os.path.splitext(__file__)[0]}"
        with open(steps_file_path, "r") as steps_file:
            steps_str = steps_file.read()
        steps = ast.literal_eval(steps_str)

        # Load the params from the file
        params_file_path = f"settings_{os.path.splitext(__file__)[0]}"
        with open(params_file_path, "r") as params_file:
            params_str = params_file.read()
        params = ast.literal_eval(params_str)

        # Run the functions in order
        for step in steps:
            method_name = step[0]
            method_params = params.get(method_name, {})
            getattr(self, method_name)(**method_params)














