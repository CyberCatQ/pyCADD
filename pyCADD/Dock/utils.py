import os
import re

from . import logger


def check_pdb(pdb: str) -> bool:
    """Check if the input is a valid PDB ID.

    Args:
        pdb (str): pdb ID

    Returns:
        bool: True if the input is a valid PDB ID, False otherwise
    """

    if re.fullmatch(r"^\d[0-9a-zA-Z]{3,}$", pdb):
        return True
    else:
        return False


def get_input_pdbid() -> str:
    """Get PDB ID from directory name or user input.

    Returns:
        str: PDB ID
    """

    pdbid = os.path.split(os.getcwd())[-1]

    if check_pdb(pdbid):
        return pdbid
    else:
        while True:
            try:
                pdbid = input("Input PDB ID:").strip().upper()
                if check_pdb(pdbid):
                    logger.info(f"PDB ID: {pdbid}")
                    return pdbid
                logger.warning(
                    f"Invalid PDB ID: {pdbid}. Please try again or press Ctrl+C to exit."
                )
            except KeyboardInterrupt:
                exit(0)
