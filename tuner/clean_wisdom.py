import glob
import shutil
import ast

files = glob.glob('../wisdom/*.wisdom')

remove_whitespace = True
remove_host = None

for f in files:
    # Backup wisdom file.
    backup = '{}.backup'.format(f)
    shutil.copyfile(f, backup)

    with open(backup, 'r') as f_in:
        with open(f, 'w') as f_out:

            for l in f_in.readlines():
                write_line = True
                l = l.strip()

                if remove_whitespace and l == '':
                    write_line = False
                else:
                    if remove_host is not None:
                        # Cast string to dictionary:
                        settings = ast.literal_eval(l)

                        if 'environment' in settings.keys() and remove_host in settings['environment']['host_name']:
                                write_line = False

                if write_line:
                    f_out.write('{}\n'.format(l))
