import matplotlib.animation as ani
from astropy.time import Time


def get_animation_writer(filename, fps=2, metadata={}, codec=None, bitrate=None, **kw):
    """
    Try to get an appropriate animation writer,
    given the filename provided.

    Parameters
    ----------

    filename : str
        The output filename string for the animation.

    fps : float
        Frames/second.

    metadata : dict
        Keywords to try to store as metadata in the animation file.

    kw : dict
        All other keywords will ignored.
    """

    # define some default metadata, which can be overwritten
    default_metadata = dict(
        artist="exoatlas",
        genre="science!",
        subject="stars and planets",
        copyright=Time.now().iso[:4],
    )
    metadata = metadata | default_metadata

    inputs = dict(fps=fps, metadata=metadata, codec=codec, bitrate=bitrate)

    # decide writer based on filename extension
    if ".mp4" in filename:
        try:
            writer = ani.writers["ffmpeg"](**inputs)
        except (RuntimeError, KeyError):
            message = """
            This computer seems unable to run ffmpeg
            to create `.mp4` movies. Try running
                `conda install ffmpeg`
            and trying again.
            """
            raise RuntimeError()
    else:
        try:
            writer = ani.writers["pillow"](**inputs)
        except (RuntimeError, KeyError):
            writer = ani.writers["imagemagick"](**inputs)
            raise RuntimeError("This computer seem unable to animate?")
    return writer
