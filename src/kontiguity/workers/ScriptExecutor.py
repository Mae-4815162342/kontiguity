from kontiguity.utils.imports import *
import logging

class ScriptExecutorScheduler(threading.Thread):
    def __init__(self, input_queue, script, output_queues = [], logger = None, **args):
        super(ScriptExecutorScheduler, self).__init__(**args)
        self._input_queue = input_queue
        self._output_queues = output_queues
        self._script = script
        self._logger = logger
        self.start()

    def run(self):
        script_executor = ScriptExecutor(self._script, logger = self._logger)
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
    def __init__(self, script, logger = None):
        self._script = script
        self._logger = logger or logging.getLogger("kontiguity")

    def launch_script(self, args_list, sbatch = False):
        str_args = [str(a) for a in args_list]

        if sbatch:
            cmd = ["sbatch", self._script] + str_args
            self._logger.info(f"Submitting (sbatch): {' '.join(cmd)}")
            subprocess.run(cmd)
            # TODO: add waiting
        else:
            cmd = ["bash", self._script] + str_args
            self._logger.info(f"Launching: {' '.join(cmd)}")
            subprocess.run(cmd)
            self._logger.info(f"Completed: {self._script} {str_args[:3]}{'...' if len(str_args) > 3 else ''}")