import fire # type: ignore

from spida.io.main import io_cli
from spida.P.main import P_cli
from spida.I.main import I_cli
from spida.__main__ import gen_run

class CLI(): 
    """
    Command Line Interface for SPIDA.
    
    Initialize the CLI with submodules for different functionalities.

    Submodules:
    - io: Handles input/output operations in SPIDA.
    - P: Preprocessing tools for SPIDA.
    - I: Inference tools for SPIDA.
    - generate_template: Generates a template for SPIDA analysis.
    """

    def __init__(self):     
        self.io = io_cli()
        self.P = P_cli()
        self.I = I_cli()
        self.generate_template = gen_run

if __name__ == "__main__":
    print("Firing")
    # fire.Fire(name="spida")
    fire.Fire(CLI)



# def main():
#     parser = argparse.ArgumentParser(prog="spida", description="SPIDA CLI")
#     subparsers = parser.add_subparsers(dest="command", required=True)

#     subparsers.add_parser(
#         "io", parents=[io_cli.get_parser()], help = "SPIDA IO tools",
#     )
    
#     subparsers.add_parser(
#         "P", parents=[P_cli.get_parser()], help = "SPIDA Preprocessing tools",
#     )

#     subparsers.add_parser(
#         "I", parents=[I_cli.get_parser()], help = "SPIDA Inference tools",
#     )

#     args = parser.parse_args()
#     func = args.func
#     kwargs = vars(args)
#     func(**kwargs)
