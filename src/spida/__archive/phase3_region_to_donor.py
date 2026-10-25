"""Archived in phase 3 of the IO redesign (2026-10-21): replaced by ``spida.io.layout.region_label``
with a mapping file (config ``REGION_LABELS`` or ``dataset_from(region_labels=...)``). The four donors
below now live in that mapping. Not imported anywhere.
"""
# flake8: noqa
# Originally in spida/utilities/sd_utils.py.
import logging

logger = logging.getLogger(__package__)


def _region_to_donor(reg_name: str) -> str:
    """
    Convert a region name to a donor name.

    Parameters:
    reg_name (str): Name of the region.

    Returns:
    str: Donor name derived from the region name.
    """
    if "2424" in reg_name:
        return "UCI2424"
    elif "5224" in reg_name:
        return "UCI5224"
    elif "4723" in reg_name:
        return "UCI4723"
    elif "7648" in reg_name:
        return "UWA7648"
    else:
        logger.warning(
            f"Region name {reg_name} does not match any known donor names. Returning the region name as donor name."
        )
        return reg_name
