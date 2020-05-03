import traceback, sys
from PySide.QtGui import *
from PySide.QtCore import *

class HaploSignal(QObject):
    def __init__(self):
        super(HaploSignal, self).__init__()

    finished = Signal()  # QtCore.Signal
    error = Signal(tuple)
    result = Signal(object)
    database = Signal(object)
    updatedatabases = Signal()
    deleted = Signal()
    deletedSeq = Signal()
    addedSeq = Signal()