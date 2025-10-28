from pyCADD.utils.tool import is_gaussian_available
if not is_gaussian_available():
    print("Gauss may not be installed correctly. "
          "Density module may not work properly.")