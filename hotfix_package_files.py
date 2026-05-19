import os

def update_line_in_file(file_path, line_number, new_content):
    """
    Updates a specific line in a file.

    :param file_path: Path to the file to be updated.
    :param line_number: The line number to update (1-based index).
    :param new_content: The new content to replace the line with.
    """
    if not os.path.exists(file_path):
        print(f"File {file_path} not found.")
        return

    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Check if the line number is valid
        if 1 <= line_number <= len(lines):
            lines[line_number - 1] = new_content + '\n'

            with open(file_path, 'w') as file:
                file.writelines(lines)

            print(f"Line {line_number} in {file_path} has been updated.")
        else:
            print(f"Invalid line number {line_number}. File {file_path} has {len(lines)} lines.")
    except Exception as e:
        print(f"An error occurred while updating {file_path}: {e}")

def update_lines_in_model_file(file_path):
    """
    Updates the specific import block in the model file.

    :param file_path: Path to the model.py file to be updated.
    """
    if not os.path.exists(file_path):
        print(f"File {file_path} not found.")
        return

    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        updated = False
        for i, line in enumerate(lines):
            if 'from collections import OrderedDict, MutableMapping' in line:
                # Replace the line with the corrected imports
                lines[i] = (
                    "from collections import OrderedDict\n"
                    "from collections.abc import MutableMapping\n"
                )
                updated = True
                break

        if updated:
            with open(file_path, 'w') as file:
                file.writelines(lines)
            print(f"The file {file_path} has been updated successfully.")
        else:
            print(f"No matching line found in {file_path} to update.")
    except Exception as e:
        print(f"An error occurred while updating {file_path}: {e}")

# Paths to the files
fields_file_path = 'lenv/lib/python3.12/site-packages/flask_restplus/fields.py'
fields_line_number = 17
fields_new_content = 'from werkzeug.utils import cached_property'

model_file_path = 'lenv/lib/python3.12/site-packages/flask_restplus/model.py'

api_file_path = 'lenv/lib/python3.12/site-packages/flask_restplus/api.py'
api_line_number = 24
api_new_content = 'from werkzeug.utils import cached_property'

# Perform updates
update_line_in_file(fields_file_path, fields_line_number, fields_new_content)
update_lines_in_model_file(model_file_path)
update_line_in_file(api_file_path, api_line_number, api_new_content)
