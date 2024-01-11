from abc import ABCMeta

class Ensemble(metaclass=ABCMeta):
    
    def __init__(self, file_path):
        self.type = "ensemble"
        self.file_path = file_path
