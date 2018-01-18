#!/usr/bin/env python
import os
import subprocess
import hashlib
from mako.template import Template
import shutil
import uuid
import docker
from status import *

class Runner(object):
    def __init__(self):
        self._status = WAITING

    def set(self, task, work_dir, input):
        self.task = task
        self.work_dir = work_dir
        self.input = input
        self.cmd = self._get_cmd()
        self.hash = hashlib.md5(self.cmd).hexdigest()
        self.task_dir = os.path.join(self.work_dir, self.task.__class__.__name__, self.hash) 
        if self._check_cached():
            self._status = CACHED
        else:
            self._status = PENDING

    def execute(self):
        self._status = RUNNING
        if os.path.exists(self.task_dir):
            shutil.rmtree(self.task_dir)
        os.makedirs(self.task_dir)
        for filepath in self.input:
            os.symlink(filepath, os.path.join(self.task_dir, os.path.basename(filepath)))
        with open(os.path.join(self.task_dir, '.command.sh'), 'w') as f:
            print >> f, '#!/bin/bash'
            for line in self.cmd.split('\n'):
                print >> f, line.strip()
        with open(os.path.join(self.task_dir, '.command.run'), 'w') as f:
            print >> f, '#!/bin/bash'
            print >> f, 'touch .start'
            print >> f, '/bin/bash -eu .command.sh 1> .command.out 2> .command.log'
            print >> f, 'echo $? > .exitcode'
        if self.task.docker != None:
            client = docker.from_env()
            client.containers.run(self.task.docker, '/bin/bash .command.run', detach=True,
                                  working_dir=self.task_dir, 
                                  volumes={self.work_dir: {'bind': self.work_dir, 'mode': 'rw'}})
        else:
            subprocess.Popen(['/bin/bash', '.command.run'],
                             cwd=self.task_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return

    def get_output(self):
        output = []
        if self.task.output() == []:
            for filename in os.listdir(self.task_dir):
                if filename.startswith('.'):
                    continue
                if not os.path.islink(os.path.join(self.task_dir, filename)):
                    output.append(os.path.join(self.task_dir, filename))
        else:
            for filename in self.task.output():
                output.append(os.path.join(self.task_dir, filename))
        return output

    def _get_cmd(self):
        current_dir = os.getcwd()
        dirname = uuid.uuid4().hex
        tmp_dir = os.path.join(self.work_dir, 'tmp', dirname)
        os.makedirs(tmp_dir)
        for filepath in self.input:
            if os.path.exists(os.path.join(tmp_dir, os.path.basename(filepath))):
                print 'error, same file name'
                exit()
            os.symlink(filepath, os.path.join(tmp_dir, os.path.basename(filepath)))
        os.chdir(tmp_dir)
        if self.task.run() is None:
            cmd = ''
        else:
            cmd_template = Template(self.task.run())
            cmd = cmd_template.render(self_=self.task)
        os.chdir(current_dir)
        shutil.rmtree(tmp_dir)
        return cmd

    def close(self):
        self._status = DONE

    @property
    def status(self):
        if self._status == RUNNING:
            if os.path.exists(os.path.join(self.task_dir, '.exitcode')):               
                with open(os.path.join(self.task_dir, '.exitcode')) as f:
                    exitcode = int(f.read())
                if exitcode == 0:
                    self._status = FINISHED
                else:
                    self._status = FAILED
        elif self._status == PENDING:
            if self._check_cached():
                self._status = CACHED
        return self._status

    def _check_cached(self):
        if not os.path.exists(self.task_dir):
            return False
        if not os.path.exists(os.path.join(self.task_dir, '.exitcode')):
            return False
        cached = True
        with open(os.path.join(self.task_dir, '.exitcode')) as f:
            exitcode = int(f.read())
            if exitcode != 0:
                cached = False
        start_date = os.path.getmtime(os.path.join(self.task_dir, '.start'))
        for filepath in self.input:
            input_date = os.path.getmtime(filepath)
            if start_date < input_date:
                cached = False
        return cached
