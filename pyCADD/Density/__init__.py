import os
if not os.path.exists(os.popen('which g16').read().strip()):
    print('Gaussian 16 is not installed or not in PATH.')
    exit()