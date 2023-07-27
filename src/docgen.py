import os, re, ast, sys


class MethodVisitor(ast.NodeVisitor):
    def __init__(self):
        self.current_class = None
        self.names = []
        self.class_names = []
        self.docstrings = []
        self.parameters = {}

    def visit_ClassDef(self, node):
        self.current_class = node.name
        self.generic_visit(node)
        self.current_class = None

    def visit_FunctionDef(self, node):
        method_name = node.name
        # If the method is part of a class, store the class name
        self.names.append(method_name)
        self.class_names.append(self.current_class)
        # Extract and store the docstring if present
        docstring = ast.get_docstring(node)
        self.docstrings.append(docstring)
        # Get the params and default values
        for arg in node.args.args:
            param_name = arg.arg
            param_default = None
            if arg.default is not ast.arg:
                param_default = ast.get_source_segment(arg.default)
            self.parameters[param_name] = param_default

        self.generic_visit(node)

    def visit_FunctionDef2(self, node):
        for arg in node.args.args:
            param_name = arg.arg
            param_default = None
            if arg.default is not ast.arg:
                param_default = ast.get_source_segment(arg.default)
            self.parameters[param_name] = param_default
        self.generic_visit(node)

class Method:
    def __init__(self, name, origin_class, doc):
        self.name = name
        self.origin_class = origin_class 
        self.doc = doc
        self.description = self.get_description()

    def get_description(self):
        if not self.doc: return ""

        # Reg-ex pattern to get description.
        pattern = r'^(.*?)(?:Parameters:|Parameters|Params|Params:|Returns|Attributes|$)'
        match = re.search(pattern, self.doc, re.DOTALL)
        if not match: return ""
        else: return match.group(1).strip()

    def __str__(self):
        comes_from = f"comes from {self.origin_class if self.origin_class else 'no class'}"
        doc = f"{'has' if self.doc else 'does not have'} documentation."
        return f"Method {self.name} {comes_from} and {doc}"


def scan_methods_in_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    # Parse the Python source code to get the AST
    try: 
        tree = ast.parse(content)
        print(f"Reading: {file_path}")
    except SyntaxError:
        print(f"Could not read: {file_path} because it contains a syntax error.")
        return None
    
    # Create the visitor and visit the AST
    visitor = MethodVisitor()
    visitor.visit(tree)
    # Make and return a list of methods
    methods = []
    for index in range(len(visitor.names)):
        # Remove any methods that start and end with __, 
        # except for init, which needs to be stripped to make markdown happy
        if visitor.names[index].startswith('__') and visitor.names[index].endswith('__'):
            if visitor.names[index] == "__init__":
                methods.append(Method(f"init {visitor.class_names[index]}", visitor.class_names[index], visitor.docstrings[index]))
            continue
        else:
            methods.append(Method(visitor.names[index], visitor.class_names[index], visitor.docstrings[index]))
    return methods

def create_param_table(method, file):
    # First, get the text in the doc with the params
    param_pattern = r'Parameters:(.*?)(?:Returns|Attributes|$)'
    param_pattern = r'(?i)(?:Parameters:|Parameters|Params|Params:)(.*?)(?:Returns|Attributes|$)'
    # param_pattern = r'Parameters:(.*?)(?:Returns|Attributes|\'\'\'|\"\"\")'
    params_text = re.search(param_pattern, method.doc, re.DOTALL)
    if not params_text:
        return ""
    params_text = params_text.group(1).strip()
    
    # Use regular expressions to find lines in the form "NAME : TYPE = DEFAULT" or "NAME : TYPE"
    param_pattern = r'^\s*(\w+)\s*:\s*(.*?)(?:\s*=\s*(.*))?$'
    # Using re.findall() to extract the name and types from each line
    params = re.findall(param_pattern,  params_text, re.MULTILINE)
    if not params:
        return ""

    # To find the description, get each line in the params text
    lines = params_text.split('\n')
    # Construct the Markdown table header
    file.write(f'\n### Parameters:\n')
    table_header = "| Name | Type | Description | Default |\n"
    table_header += "| --- | --- | --- | --- |\n"
    # Construct the table rows
    table_rows = ""
    for param in params:
        # Extract info from regex search
        name, data_type, default_value = param
        # Get the description by finding the next line.
        description = ""
        for index, line in enumerate(lines):
            if line.strip().startswith(name):
                # Check to see if the next lines are indented more than this line
                prev_spaces = len(line) - len(line.lstrip())
                next_index = index + 1
                # Loop through the bext lines
                while next_index < len(lines):
                    next_spaces = len(lines[next_index]) - len(lines[next_index].lstrip())
                    # Add the description line
                    if next_spaces > prev_spaces:
                        if description == "":
                            description = lines[next_index]
                        else: 
                            description = description + '<br />' + lines[next_index]
                        next_index = next_index + 1
                    # Break since the indentation is greater
                    else:
                        break
                # We have found (or tried to) the description of the param, so break
                break
                
        # remove any | from the data_type string
        data_type = data_type.replace(' | ', ', ')
        table_rows += f"| {name} | {data_type} | {description if description != '' else 'None'} | {default_value if default_value else 'required'} |\n"
    # Combine the header and rows to form the complete Markdown table
    markdown_table = f"\n\n{table_header}{table_rows}"
    file.write(markdown_table)
        
def generate_markdown_file(file_path, methods):
    # Create the markdown file
    base_name = os.path.basename(file_path)
    markdown_path = os.path.join('documentation', f'{os.path.splitext(base_name)[0]}.md')
    # Open the markdown file for writing
    with open(markdown_path, 'w') as markdown_file:
        # Write table of contents
        markdown_file.write('# Table of Contents\n')
        # Keep track of the class the method belongs to check if you need to indent
        current_class = None
        for method in methods:
            # Print the class name
            if method.origin_class != current_class:
                markdown_file.write(f'- {method.origin_class} Class\n')
                current_class = method.origin_class
            # Print the link to the method documentation
            prefix = " - " if not method.origin_class else "    - "
            markdown_file.write(f'{prefix}[{method.name}](#{method.name.lower().replace(" ", "-")}){": " if method.description != "" else "" }{method.description}\n')
        
        # Write documentation for each method
        for method in methods:
            # Write the name of the method
            markdown_file.write(f'\n## {method.name}\n')
            # Write info from the method documentation
            if method.doc: 
                # Write the description of the method.
                if method.description != "":
                    markdown_file.write(f'\n### Description:\n')
                    markdown_file.write(f'{method.description}\n\n')
                # Make a table for the method params
                create_param_table(method, markdown_file)

def make_documentation(directory_path):
    # Create the "documentation" directory if it doesn't exist
    os.makedirs('documentation', exist_ok=True)
    # Iterate through all Python files in the directory and the subdirectories
    for root, dirs, files in os.walk(directory_path):
        for file_name in files:
            if file_name.endswith('.py'):
                file_path = os.path.join(root, file_name)
                # Get the methods in each file
                methods = scan_methods_in_file(file_path)
                if methods:
                    # Sort the list first by the origin_class attribute with None values coming first,
                    # then alphabetically by the name attribute     
                    methods = sorted(methods, key=lambda m: (m.origin_class is not None, m.origin_class, m.name))
                    generate_markdown_file(file_path, methods)


num_args = len(sys.argv)
if num_args == 1:
    print("Invalid Input.")
    print("Please indicate the name of the directory to scan.")
else:
    directory_path = sys.argv[1]
    make_documentation(directory_path)

 

