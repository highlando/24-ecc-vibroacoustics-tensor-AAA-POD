import time


class Timer(object):
    def __init__(self, name=None, logger=None, timerinfo={}, verbose=True):
        self.name = name
        self.logger = logger
        self.timerinfo = timerinfo
        self.verbose = verbose
        self.exitmessage = ''

    def __enter__(self):
        self.tstart = time.time()
        # self.ptstart = time.process_time()

    def __exit__(self, type, value, traceback):
        elt = time.time() - self.tstart
        # elpt = time.process_time() - self.ptstart
        self.timerinfo.update(dict(elt=elt))
        if self.logger is not None:
            self.logger.info('{0}: Elapsed time: {1}'.
                             format(self.name, elt))
        if self.verbose:
            print('{0}: Elapsed time: {1}'.format(self.name, elt))
        else:
            self.exitmessage = '{0}: Elapsed time: {1}'.format(self.name, elt)
