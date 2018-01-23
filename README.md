# nextPyflow
nextPyflow is a Python module that helps you build complex bioinformatics workflow.

**Note: This project is still in alpha stage.**

## Task structure
```
class taskA(nextPyflow.Task):
    def parameter(self, args):
        pass
    def requires(self):
        return taskB()
    def run(self):
        return '''                                                                                                             
        do something                                                                          
        '''
```

Tasks are defined by Python classes which inherit from `nextPyflow.Task` class. In a task, three methods `parameter()`, `requires()` and `run()` are usually specified. The `parameter()` method defines various parameters needed for a task to execute. Also if a task need to take arguments, you can specify as arguments of `parameter()` method instead of those of `__init__()` method. The Arguments are automatically set as instance variables. The `requires()` method defines dependencies on other tasks. The `run()` method defines the command line to execute. In the `run()` method, you can use a templates library Mako's syntax and API, and you can access to the class instance as `self_`.

## Workflow Structure
### Linear
![linear](https://github.com/oliverSI/nextPyflow/blob/master/examples/linear.png)
```
class task1(nextPyflow.Task):
    pass

class task2(nextPyflow.Task):
    def requires(self):
        return task1()

class task3(nextPyflow.Task):
    def requires(self):
        return task2()
```

### Branch-Merge
![branch_merge](https://github.com/oliverSI/nextPyflow/blob/master/examples/branch_merge.png)
```
class task1(nextPyflow.Task):
    pass

class task2(nextPyflow.Task):
    def requires(self):
        return task1()

class task3(nextPyflow.Task):
    def requires(self):
        return task1()

class task4(nextPyflow.Task):
    def requires(self):
        return task2(), task3()
```

### Scatter-Gather
![branch_merge](https://github.com/oliverSI/nextPyflow/blob/master/examples/scatter_gather.png)
```
class task1(nextPyflow.Task):
    pass

class task2(nextPyflow.Task):
    def parameter(self, i):
        pass
    def requires(self):
        return task1()

class task3(nextPyflow.Task):
    def requires(self):
        for i in range(0, 5):
            yield task2(i)
```

## Example
Example pipeline for mapping sequence files against a reference genome.
```
import nextPyflow
from nextPyflow.utils import localfile

class bwa_index(nextPyflow.Task):
    def parameter(self, ref_fasta):
        self.docker = 'kathrinklee/bwa:latest'
    def requires(self):
        return localfile(self.ref_fasta)
    def run(self):
        return '''                                                                                                                  
        <%! import os %>                                                                                                            
        bwa index ${os.path.basename(self_.ref_fasta)}                                                                              
        '''
        
class bwa(nextPyflow.Task):
    def parameter(self, ref_fasta, fastq1, fastq2):
        self.docker = 'kathrinklee/bwa:latest'
    def requires(self):
        return (localfile(self.ref_fasta),
               localfile(self.fastq1),
               localfile(self.fastq2),
               bwa_index(self.ref_fasta))
    def run(self):
        return '''                                                                                                                  
        <%! import os %>                                                                                                            
        bwa mem -M \                                                                                                                
        ${os.path.basename(self_.ref_fasta)} \                                                                                      
        ${os.path.basename(self_.fastq1)} \                                                                                         
        ${os.path.basename(self_.fastq2)} \                                                                                         
        > aln-pe.sam                                                                                                                
        '''
        
class bam_sort(nextPyflow.Task):
    def parameter(self, ref_fasta, fastq1, fastq2):
        self.docker = 'comics/samtools:latest'
    def requires(self):
        return bwa(self.ref_fasta, self.fastq1, self.fastq2)
    def run(self):
        return '''                                                                                                                  
        samtools sort \                                                                                                             
        -o aln-pe-sorted.bam \                                                                                                      
        -O BAM \                                                                                                                    
        aln-pe.sam                                                                                                                  
        samtools index aln-pe-sorted.bam                                                                                            
        '''
```
Pipelines should be run like this
```
workflow = nextPyflow.Workflow(max_core=1)
workflow.run(bam_sort('/path/to/ref_fasta', '/path/to/fastq1', '/path/to/fastq2'))
```

The command line is executed in the task directory, where the symbolic link of the output from required tasks is created. In this example, a "docker" attribute is specified in the `parameter()` method, so the tasks are executed in docker container. `localfile()` is pre-defined task for importing a local file to workflow.


