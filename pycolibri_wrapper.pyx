cimport pycolibri_classes

cdef class IndexedPatternModel:
    cdef pycolibri_classes.IndexedPatternModel *thisptr
    def __cinit__(self, str filename):
        self.thisptr = new pycolibri_classes.IndexedPatternModel(filename, False)

    def __dealloc__(self):
        del self.thisptr






