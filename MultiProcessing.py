import multiprocessing
import time
import os
import datetime


class MultiProcessing:
    def __init__(self, input_file_path, output_file_path, f):
        start_time = time.time()
        self.num_cores = multiprocessing.cpu_count()
        self.infile = open(input_file_path, 'r')
        self.outfile = open(output_file_path, 'w')
        self.f = f
        self.inq = multiprocessing.Queue()
        self.outq = multiprocessing.Queue()

        self.pin = multiprocessing.Process(target=self.input_csv, args=())
        self.pout = multiprocessing.Process(target=self.output_csv, args=())
        self.ps = [multiprocessing.Process(target=self.sum_row, args=()) for i in range(self.num_cores)]

        self.pin.start()
        self.pout.start()
        for p in self.ps:
            p.start()

        self.pin.join()
        for p in self.ps:
            p.join()
        self.pout.join()

        self.infile.close()
        self.outfile.close()
        print "Time to execute: %s" % datetime.timedelta(seconds=time.time() - start_time)

    def input_csv(self):
        for i, line in enumerate(self.infile):
            self.inq.put((i, line))
        for i in range(self.num_cores):
            self.inq.put("STOP")

    def sum_row(self):
        for i, row in iter(self.inq.get, "STOP"):
            self.outq.put((i, self.f(row)))
        self.outq.put("STOP")

    def output_csv(self):
        cur = 0
        buffer = {}
        for works in range(self.num_cores):
            for i, row in iter(self.outq.get, "STOP"):
                if i != cur:
                    buffer[i] = row
                else:
                    self.outfile.write(str(row) + os.linesep)
                    cur += 1
                    while cur in buffer:
                        self.outfile.write(str(buffer[cur]) + os.linesep)
                        del buffer[cur]
                        cur += 1


def conventional(input_file_path, output_file_path, f):
    with open(input_file_path, 'r') as infile:
        with open(output_file_path, 'w') as outfile:
            for line in infile:
                outfile.write(f(line))