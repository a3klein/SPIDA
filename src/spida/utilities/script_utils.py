# Author: Amit Klein, a3klein@ucsd.edu
# some functions imported from fishtank. 
import os
import ast
import math
import argparse
import click
from pathlib import Path

import numpy as np

def parse_dict(arg: str) -> dict:
    """Parse dictionary input in the form of key1=val1,key2=val2"""
    if arg is None:
        return {}
    return {
        key: ast.literal_eval(value)
        for key, value in (item.split("=") for item in arg.split(","))
    }


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
        getattr(namespace, self.dest)["key"] = (
            "value"  # to make sure that it is never empty!
        )

        for value in values:
            key, val = value.split("=")
            # Try to convert to appropriate types
            if val.lower() in ("true", "false"):
                getattr(namespace, self.dest)[key] = val.lower() == "true"
            elif val.replace(".", "").replace("-", "").isdigit():
                if "." in val:
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
    # TODO: depreceate eventually into all click kargs parsing
    kwargs_dict = {}
    if kwargs_list:
        for kwarg in kwargs_list:
            if "=" not in kwarg:
                raise ValueError(f"Invalid kwarg format: {kwarg}. Expected 'key=value'")
            key, value = kwarg.split("=", 1)

            # Try to convert to appropriate types
            if value.lower() in ("true", "false"):
                kwargs_dict[key] = value.lower() == "true"
            elif value.replace(".", "").replace("-", "").isdigit():
                if "." in value:
                    kwargs_dict[key] = float(value)
                else:
                    kwargs_dict[key] = int(value)
            else:
                kwargs_dict[key] = value

    return kwargs_dict

def parse_click_kwargs(item_list): 
    """Parse click-style extra args into a kwargs dict.

    Accepts either:
      --key value --other_key value
    or:
      --key=value --other=value

    Values are interpreted with ast.literal_eval when possible, with these fallbacks:
      - case-insensitive 'true'/'false' -> bool
      - comma-separated -> list of strings (if literal_eval fails)
      - otherwise -> raw string
    """
    if not item_list:
        return {}

    kwargs_dict = {}
    i = 0
    n = len(item_list)

    while i < n:
        token = item_list[i]
        if not token.startswith("--"):
            raise ValueError(f"Expected key starting with '--', got {token}")

        # support --key=value and --key value
        if "=" in token:
            key, raw_val = token[2:].split("=", 1)
            i += 1
        else:
            key = token[2:]
            if i + 1 >= n:
                raise ValueError(f"Missing value for key {token}")
            raw_val = item_list[i + 1]
            i += 2

        # Normalize boolean-like lower-case tokens
        if isinstance(raw_val, str) and raw_val.lower() in ("true", "false"):
            val = raw_val.lower() == "true"
        else:
            # Try to interpret with literal_eval (numbers, lists, dicts, quoted strings)
            try:
                val = ast.literal_eval(raw_val)
            except Exception:
                # Fallback: comma-separated -> list of strings
                if isinstance(raw_val, str) and "," in raw_val:
                    val = [v.strip() for v in raw_val.split(",")]
                else:
                    # keep raw string
                    val = raw_val

        kwargs_dict[key] = val

    return kwargs_dict

def parse_json_dict(arg: dict | Path | str) -> dict:
    import json 
    # if already a dict (when called programmatically)
    if isinstance(arg, dict):
        return arg
    elif isinstance(arg, Path):
        with open(arg, 'r') as f:
            return json.load(f)
    elif isinstance(arg, str):
        # try to parse as JSON string first
        try:
            parsed = json.loads(arg)
            if not isinstance(parsed, dict):
                raise ValueError("Parsed JSON is not a dict")
            return parsed
        except Exception:
            # not a JSON string -> treat as path to a JSON file
            if os.path.exists(arg):
                with open(arg, 'r') as f:
                    return json.load(f)
            else:
                raise click.BadParameter(
                   f"{arg} must be a dict, a JSON string, or a path to a JSON file"
                )
    else:
        raise click.BadParameter(
            f"{arg} must be a dict, a JSON string, or a path to a JSON file"
        )

class JSONParam(click.ParamType):
    """Click ParamType that accepts a dict, a JSON string, or a path to a JSON file and returns a dict."""
    name = "json"

    def convert(self, value, param, ctx):
        try:
            return parse_json_dict(value)
        except click.BadParameter as e:
            self.fail(str(e), param, ctx)
