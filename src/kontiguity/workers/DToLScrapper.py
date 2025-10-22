from kontiguity.utils.imports import *

class DToLScrapperScheduler(threading.Thread):
    def __init__(self, input_queue, output_queue, **args):
        super(DToLScrapperScheduler, self).__init__(**args)
        self._input_queue = input_queue
        self._output_queue = output_queue
        self.start()

    def run(self):
        while True:
            try:
                value = self._input_queue.get(timeout=10)
            except Empty:
                break
            if value == "DONE":
                break
            limit, offset = value
            scrapper = DToLScrapper(limit, offset)
            results = scrapper.scrap()

            for result in results:
                self._output_queue.put(result)
            
class DToLScrapper():
    def __init__(self, limit, offset):
        self.URL = "https://portal.darwintreeoflife.org/api/data_portal"
        self.URL_args = {
            "limit" : limit,
            "offset" : offset,
            "sort" : "currentStatus:asc",
            "current_class" : "kingdom"
        }

    def get_url(self):
        url = self.URL + "?"
        for arg in self.URL_args.keys():
            url += f"{arg}={self.URL_args[arg]}&"
        return url[:-1]
    
    def scrap(self):
        current_url = self.get_url()
        try:
            result = requests.get(current_url)
        except:
            print("Table not received")
        return json.loads(result.text)["results"]