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


























