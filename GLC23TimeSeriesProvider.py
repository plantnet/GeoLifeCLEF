# Author: Benjamin Deneu <benjamin.deneu@inria.fr>
#         Theo Larcher <theo.larcher@inria.fr>
#
# License: GPLv3
#
# Python version: 3.10.6

from abc import abstractmethod

class TimeSeriesProvider(object):
    def __init__(self) -> None:
        pass
        
    @abstractmethod
    def __getitem__(self, item):
        pass
    
    def __repr__(self):
        return self.__str__()
    
    @abstractmethod
    def __str__(self):
        pass
    
    def __len__(self):
        pass


if __name__ == "__main__":
    p = TimeSeriesProvider()