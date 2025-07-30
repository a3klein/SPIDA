import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

import skimage as ski
import cv2


def read_info(path: str | Path) -> pd.DataFrame:
    """
    Read the MERSCOPE ULTRA data organization file

    Parameters
    ----------
    path : str or Path
        Path to the data organization file (CSV format)

    Returns
    -------
    pd.DataFrame with columns channel, prefix, round, color, frame, zPos, file_pattern
    """

    # Read in data and get used columns
    df_data_org = pd.read_csv(path)
    df_info = df_data_org[
        ["channelName", "imageType", "imagingRound", "color", "frame", "zPos"]
    ]
    df_info = df_info.rename(
        columns={
            "channelName": "channel",
            "imagingRound": "round",
            "imageType": "prefix",
        }
    )
    # assign types
    df_info = df_info.astype(
        {
            "channel": str,
            "prefix": str,
            "round": str,
            "color": int,
            "frame": "object",
            "zPos": "object",
        }
    )
    # format round
    df_info["round"] = "_" + df_info["round"]
    df_info.loc[df_info["round"] == "_-1", "round"] = ""
    # formate frame
    df_info["frame"] = df_info["frame"].apply(
        lambda x: np.asarray(x.strip("[]").split(","), dtype=int)
    )
    df_info["zPos"] = df_info["zPos"].apply(
        lambda x: np.asarray(x.strip("[]").split(","), dtype=float)
    )
    # get file pattern
    df_info["file_pattern"] = (
        df_info["prefix"] + df_info["round"] + "_{fov:04d}.jp2.zip"
    )

    return df_info


def read_zip(
    path: str | Path,
    frames: int | list = None,
    shape: tuple = (2960, 2960),
    dtype: str = "uint8",
) -> tuple[np.ndarray, dict]:
    """
    Read MERSCOPE ULTRA archived compressed .jp2.zip file

    Parameters
    ----------
    path : str | Path
        Path to the .jp2.zip file
    frames : int | list, optional
        Frames to read from the file, by default None (read all frames)
    shape : tuple, optional
        Shape of the image, by default (2960, 2960)
    dtype : str, optional
        Data type of the image, by default "uint16"

    Returns
    -------
    np.ndarray
        Image array with shape (n_frames, height, width)
    """

    zf = zipfile.ZipFile(path)
    zip_info = zf.NameToInfo
    zip_info_list = list(zip_info.items())

    if frames is None:
        frames = np.arange(len(zip_info_list) - 1)

    config = {}
    with zf.open(zip_info_list[-1][0], "r") as stream:
        for line in stream.readlines():
            key, value = line.split(b"=")
            config[key.strip().decode()] = value.strip().decode()

    _file_pattern = Path(zip_info_list[-1][0]).with_suffix("")
    _file_pattern = str(_file_pattern) + "_{fov:04d}.jp2"

    img_stack = []
    for _f in frames:
        zip_path = _file_pattern.format(fov=_f)
        zip_img = zf.open(zip_path)
        cv_img = cv2.imdecode(
            np.frombuffer(zip_img.read(), dtype=dtype), cv2.IMREAD_UNCHANGED
        )
        img_stack.append(cv_img)
    return np.asarray(img_stack), config


def read_img(
    path: str | Path,
    frames: int | list = None,
    z_project: bool = False,
    plugin: str = None,
    **plugin_args,
) -> tuple[np.ndarray, dict]:
    """
    Read image file
    """

    path = Path(path)
    suffix = path.suffix.lower()

    if frames is not None:
        if isinstance(frames, int):
            frames = [frames]
        frames = np.array(frames)

    ### WHAT IS COLOR ORDER AND WHY DO WE NEED TO SPECIFY COLORS HERE???
    # if colors is not None:
    #     if color_order is None:
    #         # TODO: determine color order from data organization file?
    #         raise ValueError("If colors is specified, color_order must also be specified")
    #     if isinstance(colors, int) or isinstance(colors, str):
    #         colors = [colors]
    #     colors = np.array(colors).astype(color_order.dtype)

    # get frames:

    if suffix == ".zip":
        img, attrs = read_zip(path, frames=frames)
    else:
        img = ski.io.imread(path, plugin=plugin, **plugin_args)

    if z_project:
        img = img.max(axis=-3)

    attrs["stage_position"] = np.asarray(
        attrs["position"].strip("()").split(","), dtype=np.float32
    )
    return img.squeeze(), attrs


def read_fov(
    path: str | Path,
    fov: int | str,
    channels: pd.DataFrame,
    z_slices: list = None,
    z_project: bool = False,
    ref_channel: str = None,
) -> tuple[np.ndarray, dict]:
    """
    Read FOV from MERSCOPE ULTRA experiment

    Parameters
    ----------
    path : str or Path
        Path to the MERSCOPE ULTRA data directory
    fov : int or str
        FOV number to read, e.g. 102 or "102"
    channels : pd.DataFrame
        DataFrame with columns ['channel', 'prefix', 'round', 'color', 'frame', 'zPos', 'file_pattern']
        containing the channel information for the experiment
    z_slices : list, optional
        List of z-slices to read, by default None (read all z-slices)
    z_project : bool, optional
        Whether to perform a z-projection (max projection) on the z-slices, by default False

    Returns
    -------
    np.ndarray
        Image array with shape (n_channels, n_frames, height, width)
        or (n_frames, height, width) if only one channel is specified
    """

    imgs = []
    attrs = []

    if isinstance(z_slices, int):
        z_slices = [z_slices]

    path = Path(path)
    # getting the frames to read
    for idx, row in channels.iterrows():
        frames = None
        if z_slices and frames is None:
            zPos_to_frames = {z: f for z, f in zip(row["zPos"], row["frame"])}
            frames = [zPos_to_frames[z] for z in z_slices]
        if z_slices is None:
            frames = row["frame"]

        file_pattern = row["file_pattern"].format(fov=fov)
        img, attr = read_img(
            path / file_pattern,
            frames=frames,
            z_project=z_project,
        )
        imgs.append(img)
        attrs.append(attr)

    if len(imgs) > 1:
        imgs = np.stack(imgs, axis=0)
    else:
        imgs = imgs[0]

    if ref_channel is not None:
        attrs = attrs[np.where(channels["channel"] == ref_channel)[0][0]]
    else:
        attrs = attrs[0]

    ## Calculating some attributes:
    z_offsets = dict()
    for i, k in enumerate(channels.loc[0, "zPos"]):
        z_offsets[i] = k
    attrs["z_offsets"] = z_offsets

    return imgs, attrs
