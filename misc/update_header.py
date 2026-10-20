import glob
import os
import re

new_year = 2024

# These names always get their end year ++'ed.
active_developers = [
        'Chiel van Heerwaarden',
        'Thijs Heus',
        'Bart van Stratum',
        'Menno Veerman',
        'Mirjam Tijhuis',
        'Steven van der Linden']

files = glob.glob('../main/*')
files += glob.glob('../src/*')
files += glob.glob('../include/*')
files += glob.glob('../python/*')
files += ['../CMakeLists.txt']

for filename in files:

    if not os.path.isdir(filename):

        with open(filename, 'r') as f:
            lines = f.readlines()

        # Find the starting point of the copyright lines,
        # and gather the copyright statements.
        lines_in = []
        start_line = None

        for n,line in enumerate(lines):

            if 'Copyright' in line:
                lines_in.append(line)

                if start_line is None:
                    start_line = n

        n_lines = len(lines_in)

        if n_lines > 0:

            # Delete old copyright lines.
            del(lines[start_line : start_line + n_lines])

            lines_out = []
            for line in lines_in:
                # Check if active contributer:
                match = re.search(r'\d{4}-\d{4} (.+)$', line)
                if match:
                    name = match.group(1)
                else:
                    raise Exception(f'Cant extract name from {line} in {filename}')

                if name in active_developers:
                    lines_out.append( re.sub(r'(\d{4})-(\d{4})', rf'\1-{new_year}', line) )
                else:
                    lines_out.append( line )

            # Insert the new header.
            lines[start_line:start_line] = lines_out

            # Save the output.
            with open(filename, 'w') as f:
                f.writelines(lines)

        else:
            print(f'WARNING: file {filename} does not have a copyright header, skipping...')
