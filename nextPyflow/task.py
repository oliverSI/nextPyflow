#!/usr/bin/env python
import inspect

class Task(object):
    all_tasks = []
    def __new__(cls, *args, **kwargs):
        for task in cls.all_tasks:
            if task.args == args and task.kwargs == kwargs and task.__class__.__name__ == cls.__name__:
                return task
        new_task =super(Task, cls).__new__(cls, *args, **kwargs)        
        cls.all_tasks.append(new_task)
        return new_task

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        args_names = inspect.getargspec(self.parameter).args[1:len(args)+1]
        args_list = [str(value) for name, value in zip(args_names, args) if not name.startswith('_')]
        kwargs_list = ['{}={}'.format(name, value) for name, value in kwargs.items() if not name.startswith('_')]
        args_kwargs_str = ', '.join(args_list+kwargs_list)
        self.name = '{}({})'.format(self.__class__.__name__, args_kwargs_str)
        self.core = 1
        self.docker = None
        keys = inspect.getargspec(self.parameter).args[1:]
        values = list(args)
        if inspect.getargspec(self.parameter).defaults is not None:
            values += inspect.getargspec(self.parameter).defaults
        for key, value in zip(keys, values):
            setattr(self, key, value)
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.parameter(*args, **kwargs)

    def parameter(self):
        pass

    def requires(self):
        return []

    def run(self):
        return
    
    def output(self):
        return []


