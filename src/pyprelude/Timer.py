import time

class Timer:
    def __init__(self, name='timer'):
        self.name = name
        self.startv = None
        self.endv = None
        self.t = -1
        self.tot = 0

    def start(self):
        self.startv = time.time()

    def stop(self):
        self.endv = time.time()
        self.t = self.endv - self.startv
        self.tot += self.t

    def reset(self):
        self.t = -1
        self.startv = None
        self.endv = None
        self.tot = 0

    def show(self):
        print(self.name, ':\t', 'total ', self.tot, '\tlap ', self.t)



