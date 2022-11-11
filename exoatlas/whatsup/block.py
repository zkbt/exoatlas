"""
Define an observing Block. This might be a particular
semester/trimester/quarter linked to an observing cycle,
or maybe it's just that you want to know what's up for
the next week or so.
"""

from ..imports import *


class Block(Talker):
    """
    The Block object, for generating a sequence of nights.

    :param start:
        The start of the block ()

    """

    def __init__(self, start=None, finish=None, plan=None, name="thisweek"):
        """
        Parameters
        ----------
        start : astropy.time.Time
            The start of the time range to consider.
            (defaults to now)
        finish : astropy.time.Time
            The end of the time range to consider.
            (defaults to start + 7 days)
        plan : Plan
            The plan to which this block is connected.
        name : str
            A descriptive name for this time range.
        """

        Talker.__init__(self)

        # store the connection to the plan
        self.plan = plan

        try:
            assert type(start) == Time
            assert start is not None
        except AssertionError:
            start = Time.now()

        try:
            assert type(finish) == Time
            assert finish is not None
        except AssertionError:
            finish = start + 7 * u.day

        # set the start and finish times (floored down to the start of the day), as JD
        zone = self.plan.observatory.standardzone.to("day").value
        self.start = np.floor(start.jd) + zone
        self.finish = np.ceil(finish.jd) + zone

        # set the noons (before the night = 0.0) and midnights to consider (0.5)
        self.noons = (np.arange(self.start, self.finish)) * u.day
        self.midnights = self.noons + 0.5 * u.day

        # maybe these could be gotten rid of?
        resolution = 5.0 / 60.0 / 24.0
        self.times = np.arange(self.start - 1, self.finish + 2.0, resolution)

        # n = (self.finish - self.start + 2*u.day)/resolution
        # self.start + np.arange(n)*resolution - 1*u.day

        self.name = name
