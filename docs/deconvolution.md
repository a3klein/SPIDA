# Deconvolution Setup & Usage

Configure and run deconvolution to improve segmentation quality using [deconwolf](https://deconwolf.fht.org/), *[Wernersson et al.](https://www.nature.com/articles/s41592-024-02294-7)*.

## Overview

Deconvolution reduces optical blur and improves segmentation quality. This requires:

1. Installing deconwolf on your system
2. Measuring your microscope parameters
3. Generating Point Spread Functions (PSFs) for your setup
4. Configuring SPIDA

## Prerequisites

### Install Deconwolf dependencies: 
deconwolf requires [fftw3](https://www.fftw.org/index.html), [libtiff](https://libtiff.gitlab.io/libtiff/), [gsl](https://coral.ise.lehigh.edu/jild13/2016/07/11/hello/), and (libpng)[https://libpng.org/pub/png/libpng.html] to be installed on your system. 

On some slurm systems these packages may already exist as modules (you can check using `module spider list libtiff` for example). If they do not I recommend installing them from binaries, taking note of where those binaries leave and adding them to your LD_LIBRARY_PATH variable either in .bashrc, or setting them as temporary variables anytime you need to run a spida command that uses deconwolf. 


### Install Deconwolf: 

download the binary from the [deconwolf github page](https://github.com/elgw/deconwolf/releases). Then follow their [instructions](https://elgw.github.io/deconwolf/start.html#installation) to build deconwolf. Make sure that either OpenCL or OpenMP are available on the system you are installing deconwolf in to enable GPU acceleration. 

After installation make sure both binaries (`dw` and `dw_bw`) are available

### Microscope Parameters You Need

1. **Numerical Aperture (NA)**: Property of the objective lens (e.g., 1.4)
2. **Refractive Index (ni)**: refractive index for the immersion (e.g. 1.518)
2. **XY Resolution (nm/pixel)**: From microscope specifications or calibration (e.g. 100 nm)
3. **Z-Step (nm)**: Distance between Z-slices in your imaging (e.g. 1.5um)
4. **Wavelengths (nm)**: What wavelengths you imaged at (e.g., 405, 488, 561)

## Generate Point Spread Functions (PSFs)

PSFs describe how your microscope blurs light for each wavelength.

### Install PSF Generation Tool

```bash
# PSF generation is typically included with deconwolf
which dw_bw  # Check if installed
```

### Generate PSFs

Create a directory for your PSF files:

```bash
mkdir -p ~/microscope_psfs/my_setup
cd ~/microscope_psfs/my_setup
```

Generate PSF for each wavelength:

```bash
# Example: NA=1.4, 100nm/pixel, 1.5um z-step, 405nm DAPI
dw_bw \
    --resxy 100 \
    --resz 1500 \
    --lambda 405 \
    --NA 1.4 \
    --ni 1.518 \
    --nslice 29 \
    psf_dir/405_z1500_psf.tiff 

# Repeat for other wavelengths
dw_bw --resxy 100 --resz 1500 --lambda 488 --NA 1.4 --ni 1.518 --nslice 29 psf_dir/488_z1500_psf.tiff 
```

**PSF File Naming Convention:**
```
{wavelength}_z{z_step}_psf.tif

Examples:
405_z1500_psf.tif
488_z1500_psf.tif
561_z1500_psf.tif
```

Verify PSFs were created:

```bash
ls -lah ~/microscope_psfs/my_setup/
```

## Configure Deconwolf for SPIDA

Create a deconwolf configuration file:

```bash
mkdir -p ~/config
cat > ~/config/deconwolf.ini << 'EOF'
[Paths]
dw_path = /path/to/dw
dw_psf_path = /home/user/microscope_psfs/my_setup
EOF
```

## Configure SPIDA

### 1. Update .env / .json file

Add the deconwolf config path:

```bash
# .env / .json
DECONWOLF_CONFIG=/home/user/config/deconwolf.ini
```

## Run Deconvolution

### Basic Usage

TODO

## Troubleshooting

### PSF Files Not Found

```
Error: PSF file not found: 405_z500_psf.tif
```

**Solution:**
- Check PSF directory path in `deconwolf.ini`
- Verify PSF naming matches: `{wavelength}_z{z_step}_psf.tif`
- List PSF files: `ls {dw_psf_path}/`

### Wrong Microscope Parameters

```
Error: NA mismatch between PSF and parameters
```

**Solution:**
- Verify `na`, `xy_res`, `z_step` match your microscope
- Check PSF file names match your z_step
- Regenerate PSFs if parameters change

### Deconvolution Too Slow

**Solutions:**
- Use GPU if available
- Inrease / Decrease tilesize depending on where the bottleneck is: 
    - If it is in the deconwolf runtime, decrease tilesize. 
    - If it is in the time to transfer images to GPU, increase tilesize.

### Output Images Look Bad
-- Can increase / decrease the number of padding on the z-axis used for the deconvolution (edge z stacks may not come out as good).

## Output Files

Deconvolved images are saved as:

```
{ROOT_DIR}/{experiment_name}/analysis/{region_name}/tiles/
├── 405_z{#}.decon.tif        # DAPI deconvolved
├── 488_z{#}.decon.tif        # PolyT deconvolved
```

## See Also

- [S Module (Spatial Processing)](./S_USAGE.md) - Segmentation after deconvolution
- [Configuration Guide](./CONFIGURATION.md) - Full configuration setup