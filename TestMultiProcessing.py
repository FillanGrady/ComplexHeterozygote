import multiprocessing
import multiprocessing.queues
import collections
import os
import matplotlib.pyplot as plt
import numpy as np
import time


class Data:
    def __init__(self):
        self.annotated_file_path = 'AnnotatedTest.vcf'
        self.output_file_path = 'OutputTest.vcf'
        self.count_file_path = 'CountTest.vcf'
        self.annotated_file = open(self.annotated_file_path, 'r')
        self.output_file = open(self.output_file_path, 'w')
        self.count_file = open(self.count_file_path, 'w')
        self.count_dict = collections.defaultdict(int)
        self.num_procs = multiprocessing.cpu_count()
        self.inq = multiprocessing.Queue()
        self.countq = multiprocessing.Queue()
        self.count_outq = multiprocessing.Queue()
        self.outq = multiprocessing.Queue()
        self.begin = time.time()

        self.timeq = multiprocessing.Queue()

        self.pin = multiprocessing.Process(target=self.read_lines, args=())
        self.pout = multiprocessing.Process(target=self.write_lines, args=())
        self.ps = [multiprocessing.Process(target=self.process_line, args=()) for _ in range(int(self.num_procs / 2))]
        self.pcount = [multiprocessing.Process(target=self.count_line, args=()) for _ in range(int(self.num_procs / 2))]

        self.pin.start()
        self.pout.start()
        for p in self.ps:
            p.start()
        for p in self.pcount:
            p.start()
        self.pin.join()
        for p in self.ps:
            p.join()
        for p in self.pcount:
            p.join()
        self.pout.join()
        self.count_outq.put("STOP")

        self.output_file.close()
        for dictionary in iter(self.count_outq.get, "STOP"):
            self.timeq.put((time.time(), "count_outq-"))
            if isinstance(dictionary, collections.defaultdict):
                for key, value in dictionary.items():
                    self.count_dict[key] += value
        for key, value in self.count_dict.items():
            string = str(key) + "\t" + str(value) + os.linesep
            self.count_file.write(string)
        self.count_file.close()
        self.timeq.put("STOP")

        outq_size = 0
        countq_size = 0
        inq_size = 0
        count_outq_size = 0

        times = np.linspace(0, .8, 100)
        sizes = np.zeros((np.size(times), 4))
        index = 0
        for t, action in iter(self.timeq.get, "STOP"):
            t -= self.begin
            while t > times[index]:
                sizes[index, 0] = outq_size
                sizes[index, 1] = countq_size
                sizes[index, 2] = inq_size
                sizes[index, 3] = count_outq_size
                index += 1
            if action == "outq+":
                outq_size += 1
            elif action == "outq-":
                outq_size -= 1
            elif action == "countq+":
                countq_size += 1
            elif action == "countq-":
                countq_size -= 1
            elif action == "inq+":
                inq_size += 1
            elif action == "inq-":
                inq_size -= 1
            elif action == "count_outq+":
                count_outq_size += 1
            elif action == "count_outq-":
                count_outq_size -= 1
            else:
                print(action)
                raise ValueError
        out_plot, = plt.plot(times, sizes[:, 0])
        count_plot, = plt.plot(times, sizes[:, 1])
        in_plot, = plt.plot(times, sizes[:, 2])
        count_out_plot, = plt.plot(times, sizes[:, 3])
        plt.legend([out_plot, count_plot, in_plot, count_out_plot], ["Outq", "Countq", "Inq", "Count_outq"])
        plt.show()

    def read_lines(self):
        self.annotated_file.readline()
        for line in self.annotated_file:
            self.inq.put(line)
            self.timeq.put((time.time(), "inq+"))
        for _ in range(self.num_procs):
            self.inq.put("STOP")

    def process_line(self):
        for line in iter(self.inq.get, "STOP"):
            words = line.strip().split("\t")
            self.timeq.put((time.time(), "inq-"))
            self.countq.put(words[9:])
            self.timeq.put((time.time(), "countq+"))
            self.timeq.put((time.time(), "outq+"))
            self.outq.put(words[1:9])
        self.countq.put("STOP")
        self.outq.put("STOP")

    def count_line(self):
        temp_dict = collections.defaultdict(int)
        for line in iter(self.countq.get, "STOP"):
            self.timeq.put((time.time(), "countq-"))
            for index, word in enumerate(line):
                temp_dict[index] += int(word[0])
                temp_dict[index] += int(word[-1])
        self.count_outq.put(temp_dict,)
        self.timeq.put((time.time(), "count_outq+"))

    def write_lines(self):
        for line in iter(self.outq.get, "STOP"):
            self.output_file.write("\t".join(line) + os.linesep)
            self.timeq.put((time.time(), "outq-"))


if __name__ == "__main__":
    Data()