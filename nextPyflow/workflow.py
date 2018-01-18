#!/usr/bin/env python
import os
import networkx as nx
import logging
import shutil
import time
import signal
from status import *
import runner

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('[%(levelname)s]: %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

class Workflow(object):
    def __init__(self, work_dir='work', max_core=1, params=None):
        self.work_dir = os.path.abspath(work_dir)
        self.max_core = max_core
        self.processes = []
        self.core = 0
        self.params = params
        
    def _make_dag(self, task, dag=nx.DiGraph()):
        if hasattr(task.requires(), '__iter__'):
            requires = task.requires()
        else:
            requires = [task.requires()]
        for required_task in requires:
            dag.add_node(required_task, runner=runner.Runner())
            dag.add_node(task, runner=runner.Runner())
            dag.add_edge(required_task, task)
            self._make_dag(required_task, dag)
        return dag
    
    def _get_executable_tasks(self):
        executable_tasks = []
        for task in self.dag.nodes: 
            runner = self.dag.nodes[task]['runner']
            if not(runner.status == WAITING or runner.status == PENDING):
                continue
            required_tasks = self.dag.predecessors(task)
            required_runners = [self.dag.nodes[required_task]['runner'] for required_task in required_tasks]
            statuses = [required_runner.status == DONE for required_runner in required_runners]
            if all(statuses):
                executable_tasks.append(task)
        return executable_tasks
    
    def start(self):
        for task in self._get_executable_tasks():
            runner = self.dag.nodes[task]['runner']
            if runner.status == WAITING:
                input = []
                for required_task in self.dag.predecessors(task):
                    required_runner =  self.dag.nodes[required_task]['runner']
                    input += required_runner.get_output()
                runner.set(task, self.work_dir, input)            
            if runner.status == PENDING:
                if self.max_core >= self.core + task.core:
                    logger.info('Starting %s, hash:%s', task.name, runner.hash)
                    runner.execute()
                    self.core += task.core
        for task in self.dag.nodes:
            runner = self.dag.nodes[task]['runner']
            if runner.status == FINISHED:
                runner.close()
                self.core -= task.core
                logger.info('Finished %s', task.name)
            elif runner.status == CACHED:
                runner.close()
                logger.info('Using cache %s, hash:%s', task.name, runner.hash)
            elif runner.status == FAILED:
                logger.error('%s', task.name)
                logger.error('%s', os.path.join(runner.task_dir, '.command.log'))
                exit()        
        statuses = [self.dag.nodes[task]['runner'].status == DONE for task in self.dag.nodes]        
        return all(statuses)
            
    def run(self, task):
        self.dag = self._make_dag(task)
        signal.signal(signal.SIGINT, self._exit)
        while True:
            if self.start():
                break
            time.sleep(1)
        tmp_dir = os.path.join(self.work_dir, 'tmp')
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)

    def _exit(self, *args):
        logger.error('Exiting')
        tmp_dir = os.path.join(self.work_dir, 'tmp')
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        for task in self.dag.nodes:
            runner = self.dag.nodes[task]['runner']
            if runner.status == RUNNING:
                shutil.rmtree(self.dag.nodes[task]['runner'].task_dir)
        exit()


