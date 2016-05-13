import glob

class _Empty: pass

# Convert a string to int/float/str
def _int_or_float_or_str(value):
    try:
        if ('.' in value or 'e' in value):
            return float(value)
        else:
            return int(value)
    except:
        return value.rstrip()

# Convert namelist value or list
def _convert_value(value):
    if (',' in value):
        value = value.split(',') 
        return [_int_or_float_or_str(val) for val in value]
    else:
        return _int_or_float_or_str(value)

# Read namelist
class Read_namelist:
    def __init__(self, namelist_file=None):
        if (namelist_file is None):
            namelist_file = glob.glob('*.ini')[0]

        with open(namelist_file) as f:
            for line in f:
                if (len(line.strip()) > 0):
                    if (line.strip()[0] == '[' and line.strip()[-1] == ']'):
                        curr_group_name = line.strip()[1:-1]
                        curr_group = _Empty()
                        setattr(self, curr_group_name, curr_group)
                    else:
                        setattr(curr_group, line.split('=')[0], _convert_value(line.split('=')[-1]))

# Get the cross-section indices
def get_cross_indices(variable, mode):
    files = glob.glob('{}.{}.*.*'.format(variable, mode))
    time = files[0].split('.')[-1]
    files = glob.glob('{}.{}.*.{}'.format(variable, mode, time))
    return [int(f.split('.')[-2]) for f in files]
