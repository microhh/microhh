import glob
import os

files = glob.glob('../main/*')
files += glob.glob('../src/*')
files += glob.glob('../include/*')
files += glob.glob('../python/*')
files += ['../CMakeLists.txt']

for filename in files:

    if not os.path.isdir(filename):

        fileread = open(filename, "r")
        lines = fileread.readlines()
        fileread.close()

        # Find the line number where the copyright info starts.
        found_copyright = False

        nline = 0
        for n in lines:
            # Store the commenting style, to make sure C++, Python and Cmake work.
            npos = n.find('Copyright')
            if (npos != -1):
                found_copyright = True
                left_of_copyright = n[0:npos]
                break
            nline += 1

        if found_copyright:

            # Delete the three lines of this header.
            del(lines[nline:nline + 3])

            newlines = [f'{left_of_copyright}Copyright (c) 2011-2024 Chiel van Heerwaarden\n',
                        f'{left_of_copyright}Copyright (c) 2011-2024 Thijs Heus\n',
                        f'{left_of_copyright}Copyright (c) 2014-2024 Bart van Stratum\n']

            # Insert the new header.
            lines[nline:nline] = newlines[:]

            # Save the output.
            filewrite = open(filename, "w")
            filewrite.writelines(lines)
            filewrite.close()

        else:
            print(f'WARNING: file {filename} does not have a copyright header, skipping...')
