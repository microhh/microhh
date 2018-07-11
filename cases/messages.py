#
# This file is part of LS2D.
#
# Copyright (c) 2017-2018 Bart van Stratum
#
# LS2D is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# LS2D is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with LS2D.  If not, see <http://www.gnu.org/licenses/>.
#

import datetime

_opts = {
   'blue'   : '\033[94m',
   'green'  : '\033[92m',
   'purple' : '\033[95m',
   'red'    : '\033[91m',
   'yellow' : '\033[93m',
   'bf'     : '\033[1m',
   'ul'     : '\033[4m',
   'end'    : '\033[0m'
}

def microhh():
    print(" __  __ _                _   _ _   _ ")
    print("|  \/  (_) ___ _ __ ___ | | | | | | |")
    print("| |\/| | |/ __| '__/ _ \| |_| | |_| |")
    print("| |  | | | (__| | | (_) |  _  |  _  |")
    print("|_|  |_|_|\___|_|  \___/|_| |_|_| |_|")

def print_header(message, time=True):
    """
    Format of print statements indicating new main routine
    """
    if time:
        now = datetime.datetime.now()
        print('[{}] {}{}{}'.format(now.strftime('%d-%m: %H:%M'), _opts['green'], message, _opts['end']))
    else:
        print('{}{}{}{}'.format(_opts['green'], _opts['bf'], message, _opts['end']))

def print_message(message):
    """
    Format of print statements
    """
    print(' - {}'.format(message))

def print_warning(message):
    """
    Format of print warnings
    """
    print('{}{}WARNING:{} {}'.format(_opts['yellow'], _opts['bf'], _opts['end'], message))

def print_error(message):
    """
    Format of print errors
    """
    print('{}{}ERROR:{} {}'.format(_opts['red'], _opts['bf'], _opts['end'], message))


if __name__ == '__main__':
    microhh()
    print_header('Running task ASDF')
    print_warning('This might not be okay')
    print_error('This is not okay')
