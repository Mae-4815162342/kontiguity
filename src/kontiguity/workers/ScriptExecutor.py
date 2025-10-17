from kontiguity.utils.imports import *

class ScriptExecutorScheduler(threading.Thread):
    def __init__(self, input_queue, script, output_queues = [], **args):
        super(ScriptExecutorScheduler, self).__init__(**args)
        self._input_queue = input_queue
        self._output_queues = output_queues
        self._script = script
        self.start()

    def run(self):
        script_executor = ScriptExecutor(self._script)
        while True:
            try:
                value = self._input_queue.get(timeout= 10)
            except Empty:
                break
            if value == "DONE":
                break
            args, sbatch = value
            output = script_executor.launch_script(args, sbatch=sbatch)
            for queue in self._output_queues:
                queue.put(output)

class ScriptExecutor():
    def __init__(self, script):
        self._script = script

    def launch_script(self, args_list, sbatch = False):
        if sbatch:
            subprocess.run(["sbatch", self._script] + args_list)
        else:
            subprocess.run(["bash", self._script] + args_list)