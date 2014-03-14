"""schedule.py - cMonkey scheduling module

Scheduling is a central feature in cMonkey and this module attempts to
make it simple and flexible to create schedules. Schedules are typically
defined in a .ini file.
"""


class RepeatingSchedule:
    """A basic building block of a schedule: start and interval"""

    def __init__(self, start, interval):
        self.start = start
        self.interval = interval

    def __call__(self, iteration):
        return iteration >= self.start and (iteration - self.start) % self.interval == 0

    def __repr__(self):
        return '%d,%d' % (self.start, self.interval)

    def __str__(self):
        return repr(self)

class OneTimeSchedule:
    """A basic building block of a schedule: runs only in one iteration"""

    def __init__(self, iteration):
        self.iteration = iteration

    def __call__(self, iteration):
        return self.iteration == iteration

    def __repr__(self):
        return '%d' % iteration

    def __str__(self):
        return repr(self)


class CompositeSchedule:
    """A composite of one or more schedules. It asks its sub schedules"""

    def __init__(self, schedules):
        self.schedules = schedules

    def __call__(self, iteration):
        for schedule in self.schedules:
            if schedule(iteration):
                return True
        return False

    def __repr__(self):
        return ':'.join(map(str, self.schedules))

    def __str__(self):
        return repr(self)


def make_schedule(schedulestr):
    """creates a schedule for the specified schedule string.
    The following formats is supported
    <schedule>:[<schedule>]*
    where schedule is one of

    start,interval - repeating
    iteration - one-time
    """
    def make_sched(s):
        if ',' in s:
            start, interval = tuple(map(int, s.split(',')))
            return RepeatingSchedule(start, interval)
        else:
            return OneTimeSchedule(int(s))

    scheds = schedulestr.split(':')
    if len(scheds) > 1:
        return CompositeSchedule([make_sched(sched) for sched in scheds])
    else:
        return make_sched(scheds[0])
