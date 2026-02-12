#!/usr/bin/env python3
"""
CLI entrypoint for SPIDA P module using click.
"""

import logging

from spida.settings import configure_logging_for_runtime
from spida.P.cli import cli


if __name__ == "__main__":
    # Configure root logger with INFO level handlers (to allow INFO messages through)
    env = configure_logging_for_runtime(
        level=logging.INFO,  # Handlers need to accept INFO level
    )

    # Set root logger level to WARNING to suppress other modules
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.WARNING)

    # Set the entire spida package to INFO level as well
    spida_logger = logging.getLogger("spida")
    spida_logger.setLevel(logging.INFO)

    # Configure the spida.P module logger (parent for all files in this module)
    module_logger = logging.getLogger("spida.P")
    module_logger.setLevel(logging.INFO)

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    logger.info(f"Logging configured for environment: {env}")

    cli()
