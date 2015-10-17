import heapq
import itertools

class PQueue:
    _REMOVED = '<removed-task>'
    def __init__(self):
        self._heap = []
        self._entry = {}
        self._counter = itertools.count()

    def add(self, task, priority=0):
        'Add a new task or update the priority of an existing task'
        if task in self._entry:
            self.remove(task)
        count = next(self._counter)
        entry = [priority, count, task]
        self._entry[task] = entry
        heapq.heappush(self._heap, entry)

    def remove(self, task):
        'Mark an existing task as REMOVED.  Raise KeyError if not found.'
        entry = self._entry.pop(task)
        entry[-1] = self._REMOVED

    def pop(self):
        'Remove and return the lowest priority task. Raise KeyError if empty.'
        while self._heap:
            priority, count, task = heapq.heappop(self._heap)
            if task is not self._REMOVED:
                del self._entry[task]
                return task
        raise KeyError('pop from an empty priority queue')

    def filterGet(self, filterFun):
        ret = []
        for entry in self._heap:
            if entry[-1] is not self._REMOVED and filterFun(entry[-1]):
                ret[:0] = [entry[-1]]
        return ret
