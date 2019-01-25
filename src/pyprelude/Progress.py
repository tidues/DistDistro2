import time

class Progress:
    def __init__(self, total_iter, response_time = 60):
        self.tot = total_iter * 1.0
        self.rt = response_time
        self.cnt = 0
        self.start = None

    def count(self):
        self.cnt += 1
        if self.start is None:
            self.start = time.time()
        else:
            end = time.time()
            if end - self.start > self.rt:
                print('progress:\t%0.4f%%' % (self.cnt/self.tot * 100))
                self.start = end




