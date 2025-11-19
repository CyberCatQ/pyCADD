import importlib
import subprocess
import sys
import os

from . import logger


def check_virtual_environment(package_name, environment_name):
    try:
        importlib.import_module(package_name)
        logger.debug(f"Already in virtual environment '{environment_name}'.")
        return True
    except ImportError:
        logger.debug(f"Not in virtual environment '{environment_name}'.")
        return False


def create_environment(environment_name, schrodinger_home, root_path=None):
    if root_path is None:
        root_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    env_path = os.path.join(root_path, environment_name)
    if os.path.exists(env_path):
        try:
            import pyCADD

            logger.debug(f"Virtual environment '{environment_name}' already exists.")
            return env_path
        except ImportError:
            pass

    logger.info("Initializing Schrodinger virtual environment for first runing ...")
    try:
        subprocess.run(
            f"{schrodinger_home}/run schrodinger_virtualenv.py {env_path}",
            shell=True,
            check=True,
            stdout=subprocess.DEVNULL,
        )
        subprocess.run(
            f"bash -c 'source {env_path}/bin/activate && pip install -e {root_path}'",
            shell=True,
            check=True,
            stdout=subprocess.DEVNULL,
        )
        logger.info(f"Schrodinger virtual environment created successfully at {env_path}.")
        return env_path
    except subprocess.CalledProcessError as e:
        print(f"Error creating virtual environment: {e}")
        sys.exit(1)


def _save_schrodinger_dir(schrodinger_home: str):
    if not os.path.exists(os.path.join(schrodinger_home, "run")):
        print(f"\033[31mError: {schrodinger_home} is not a valid Schrodinger installation.\033[0m")
        sys.exit(1)
    else:
        with open(os.path.join(os.environ.get("HOME"), ".schrodinger_dir"), "w") as f:
            f.write(schrodinger_home)


def _get_schrodinger_dir():
    try:
        with open(os.path.join(os.environ.get("HOME"), ".schrodinger_dir"), "r") as f:
            schrodinger_home = f.read().strip()
    except FileNotFoundError:
        schrodinger_home = None
    return schrodinger_home


def main():
    args = sys.argv[1:]
    if len(args) > 0:
        args = " ".join(args)

    package_name = "schrodinger"
    environment_name = "pycadd-dock.ve"
    schrodinger_home = os.environ.get("SCHRODINGER") or _get_schrodinger_dir()

    if not schrodinger_home:
        print("SCHRODINGER environment variable is not set.")
        schrodinger_home = input("Please enter the path to Schrodinger installation: ")
        _save_schrodinger_dir(schrodinger_home)

    if not check_virtual_environment(package_name, environment_name):
        env_path = create_environment(environment_name, schrodinger_home)
        logger.debug("Re-running in the virtual environment.")
        subprocess.run(
            f"bash -c 'source {env_path}/bin/activate > /dev/null && pycadd-dock {args}'",
            shell=True,
            check=True,
        )
        sys.exit(0)

    from pyCADD.Dock.cli import cli_main

    cli_main()


if __name__ == "__main__":
    main()
