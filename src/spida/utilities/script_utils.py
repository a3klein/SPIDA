# From Fishtank
import ast
import math
import argparse
from pathlib import Path

import numpy as np


def parse_dict(arg: str) -> dict:
    """Parse dictionary input in the form of key1=val1,key2=val2"""
    if arg is None:
        return {}
    return {key: ast.literal_eval(value) for key, value in (item.split("=") for item in arg.split(","))}


def parse_list(arg: str) -> list:
    """Parse list input separated by commas"""
    if arg is None:
        return None
    return arg.split(",")


def parse_index(arg: str) -> list:
    """Parse index input in the form of start:end:step"""
    if arg is None:
        return None
    if "," in arg:
        return parse_list(arg)
    if ":" in arg:
        return list(range(*map(int, arg.split(":"))))
    return [int(arg)]


def parse_path(arg: str) -> Path:
    """Parse path input"""
    if arg is None:
        return None
    return Path(arg)


def parse_bool(arg: str) -> bool:
    """Parse boolean input"""
    if arg is None:
        return False
    return arg.lower() in ["true", "1", "t", "y", "yes"]


def parse_rotation(arg: str) -> float:
    """Parse rotation input"""
    if arg is None:
        return None
    elif arg.endswith(".npy"):
        matrix = np.load(Path(arg))
        radians = math.atan2(matrix[1, 0], matrix[0, 0])
        return math.degrees(radians)
    return float(arg)


# Parssing kwargs in the argparse call as an action 
class ParseKwargs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        getattr(namespace, self.dest)["key"] = "value" # to make sure that it is never empty!

        for value in values:
            key, val = value.split('=')
            # Try to convert to appropriate types
            if val.lower() in ('true', 'false'):
                getattr(namespace, self.dest)[key] = val.lower() == 'true'
            elif val.replace('.', '').replace('-', '').isdigit():
                if '.' in val:
                    getattr(namespace, self.dest)[key] = float(val)
                else:
                    getattr(namespace, self.dest)[key] = int(val)
            else:
                getattr(namespace, self.dest)[key] = val

# parsing kwargs in the main blcok as a functino
def parse_kwargs(kwargs_list):
    """
    Parse a list of key=value strings into a dictionary.

    Args:
        kwargs_list: List of strings in 'key=value' format
        
    Returns:
        dict: Dictionary of parsed key-value pairs
    """
    kwargs_dict = {}
    if kwargs_list:
        for kwarg in kwargs_list:
            if '=' not in kwarg:
                raise ValueError(f"Invalid kwarg format: {kwarg}. Expected 'key=value'")
            key, value = kwarg.split('=', 1)
            
            # Try to convert to appropriate types
            if value.lower() in ('true', 'false'):
                kwargs_dict[key] = value.lower() == 'true'
            elif value.replace('.', '').replace('-', '').isdigit():
                if '.' in value:
                    kwargs_dict[key] = float(value)
                else:
                    kwargs_dict[key] = int(value)
            else:
                kwargs_dict[key] = value
    
    return kwargs_dict