# Software Setup Guide

This file will help you set up the necessary tools to be able to run through the standard SPIDA pipeline to preprocess MERSCOPE data. 

## Deconvolution 

Deconvolution is used to increase the accuracy / performance of the segmentation. It is highly recommended for brain tissue which lacks a good cell boundary signal, or for regions like the dentate gyrus / cerebellum which contain very dense cells. It might not be necessary for other tissues. 

Refer to the [Decovolution](./deconvolution.md) page for instructions on how to set up the deconvolution pipeline. 

## VPT (vizgen-Postprocessing-tool) — optional

> **You do not need VPT for normal use.** The standard post-processing steps (cell-by-gene
> matrix, metadata, image-intensity, geometry conversion) are reimplemented in pure Python
> and run by default (`--backend native`). VPT is kept **only as a backwards-compatibility
> fallback** (`--backend vpt`) for pre-redesign workflows/data, and it is **slower** than
> native. Outputs are equivalent up to small library-version differences (e.g. GDAL
> `rasterize` affects `sum_signals` by ~0.5–1%; shapely version affects `anisotropy`). You
> can skip this whole section unless you specifically need the vpt backend.

VPT reference documentation: [vizgen_postprocessing](https://vizgen.github.io/vizgen-postprocessing/index.html).

If you do want the optional vpt backend, install it using poetry into a fresh conda
environment. You can follow the vpt instructions, but here is the tldr of how I do it: 
```
micromamba create -n vpt python=3.10 poetry -y
micromamba activate vpt
cd /path/to/programs/ # path to wherever you want to put the clones directory (can be removed afterwards)
git clone https://github.com/Vizgen/vizgen-postprocessing
cd vizgen-postprocessing
poetry install --all-extras
rm -rf ~/.cache/pypoetry # optional step to remove .cache
rm -rf /path/to/programs/vizgen-postprocessing # optional step to remove the dir
```

Once you have a binary file which contains vpt make sure to update the `VPT_BIN_PATH` parameter of the configurations file (or the .env file in the base SPIDA folder) with the full path to that bin. In the above setup the binary path will be `/path-to-micromamba_dir/envs/vpt/bin`

## Proseg

[Proseg](https://github.com/dcjones/proseg), *[D.C. Jones, et al.](https://www.nature.com/articles/s41592-025-02697-0)* is a transcript aware segmentation algorithm that attempts to optimize cell segmentation seeds. In this package we offer the option to use proseg on top of Cellpose to enhance segmentations. 

Proseg runs on rust, so to install it you need to first install Rust an Cargo ([guide](https://doc.rust-lang.org/cargo/getting-started/installation.html)), and then use cargo to install proseg. I would make sure that your version of proseg is > 3.10.

```
curl https://sh.rustup.rs -sSf | sh
source ~/.bashrc
cargo install proseg
```

Your `RUST_BIN_PATH` should then be under `~/.cargo/bin`. 

## Configuring .env file

Typically, relative filepaths for projects can change, however these software dependencies stay the same. I will at the very least make sure that the .env file in the SPIDA root project folder is updated with the `DECONWOLF_CONFIG` and `RUST_BIN_PATH`, since these parameters are unlikely to change between projects on the same machine. (`VPT_BIN_PATH` is optional — only set it if you use the `--backend vpt` fallback.) 

You can refer to the [Configuration](./configuration.md) page for more information.