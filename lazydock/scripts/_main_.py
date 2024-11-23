import importlib
import os
import sys

from mbapy_lite.base import get_fmt_time
from mbapy_lite.file import opts_file

if __name__ == '__main__':
    from lazydock.utils import get_storage_path
else:
    from ..utils import get_storage_path

scripts_info = opts_file(get_storage_path('lazydock-cli-scripts-list.json'), way = 'json')
exec2script = {exec_name: script_info['script name'] for script_info in scripts_info.values() for exec_name in script_info['exec_names']}


def print_version_info():
    import lazydock
    print('lazydock python package command-line tools')
    print('lazydock version: ', lazydock.__version__, ', build: ', lazydock.__build__)
    print('lazydock author: ', lazydock.__author__, ', email: ', lazydock.__author_email__)
    print('lazydock url: ', lazydock.__url__, ', license: ', lazydock.__license__)
    
def print_help_info():
    help_info = """
            usage-1: lazydock-cli [-h] [-l | -i]
            options:
            -h, --help  show this help message and exit
            -l, --list  print scripts list
            -i, --info  print scripts info
            usage-2: lazydock-cli [sub-scripts-name] [args] [-h]
            options:
            sub-scripts-name  name of scripts in lazydock.scripts
            args  args for sub-scripts
            -h, --help  show this help message and exit
            """
    print(help_info)
    
def print_scripts_list():
    for idx, script in enumerate(scripts_info):
        print(f'scripts {idx:3d}: {scripts_info[script]["name"]}')
        print(f'script file name: {scripts_info[script]["script name"]}')
        print(f'exec names: {", ".join(scripts_info[script]["exec_names"])}')
        print(scripts_info[script]['brief'])
        print('-'*100)

def print_scripts_info():
    for idx, script in enumerate(scripts_info):
        print(f'scripts {idx:3d}: {scripts_info[script]["name"]}')
        print(f'script file name: {scripts_info[script]["script name"]}')
        print(f'exec names: {", ".join(scripts_info[script]["exec_names"])}')
        print(scripts_info[script]['brief'])
        print(scripts_info[script]['detailed'])
        print('-'*100)
        
def exec_scripts():
    import lazydock

    # check --pause-after-exec argumet
    pause_after_exec = '--pause-after-exec' in sys.argv
    if pause_after_exec:
        sys.argv.remove('--pause-after-exec')
    # check and exec scripts NOTE: DO NOT use exec
    script_name = exec2script[sys.argv[1]]
    script = importlib.import_module(f'.{script_name}', 'lazydock.scripts')
    script.main(sys.argv[2:])
    # pause if --pause-after-exec
    if pause_after_exec:
        os.system('pause') # avoid cmd window close immediately
    
def main():  
    def _handle_unkown():
        print(f'lazydock-cli: unkown scripts: {sys.argv[1]} and args: {", ".join(sys.argv[2:])}, SKIP.')
        print('bellow are all scripts list:\n\n')
        print_scripts_list()
        
    if len(sys.argv) == 1:
        print_version_info()
    elif len(sys.argv) == 2:
        if sys.argv[1] in ['-l', '--list']:
            print_scripts_list()
        elif sys.argv[1] in ['-i', '--info']:
            print_scripts_info()
        elif sys.argv[1] in ['-h', '--help']:
            print_help_info()
        elif sys.argv[1] in exec2script:
            # exec scripts with only ZERO arg
            exec_scripts()
        elif os.path.exists(sys.argv[1]) and sys.argv[1].endswith('.mpss'):
            # load lazydock Script Session file and exec
            from lazydock.scripts._script_utils_ import Command
            print(f'loading session from file: {sys.argv[1]}')
            Command(None).exec_from_session(sys.argv[1])
            os.system('pause') # avoid cmd window close immediately
        else:
            _handle_unkown()
    else:
        if sys.argv[1] in exec2script:
            # exec scripts
            exec_scripts()
        else:
            _handle_unkown()
    # exit
    print(f'\nlazydock-cli: exit at {get_fmt_time()}')
            

if __name__ == '__main__':
    main()