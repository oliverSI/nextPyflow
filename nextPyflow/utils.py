#!/usr/bin/env python
import task

class localfile(task.Task):
    def parameter(self, filename):
        self.filename = filename

    def run(self):
        return '''
        cp -rL ${self_.filename} .
        '''



