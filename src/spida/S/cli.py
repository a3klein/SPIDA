import argparse 
import spida.S.decon_script as decon_script

def main(): 
    """ Main Entry Point for the S module os SPIDA"""

    parser = argparse.ArgumentParser(prog="spida.S", description="SPIDA S module CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: decon_image
    decon_parser = subparsers.add_parser(
        'decon_image', parents=[decon_script.get_parser()], help='Deconvolve large image files in tiles'
    )

    args = parser.parse_args()
    func = args.func
    kwargs = vars(args)
    func(**kwargs)

if __name__ == "__main__":
    main()