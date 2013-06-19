cimport pycolibri_classes

cdef class IndexedPatternModel:
    cdef pycolibri_classes.IndexedPatternModel *thisptr
    def __cinit__(self, str filename):
          self.thisptr = new pycolibri_classes.IndexedPatternModel(filename, False)
    def __dealloc__(self):
        del self.thisptr


cdef class ClassEncoder:
    cdef pycolibri_classes.ClassEncoder *thisptr
    def __cinit__(self, str filename):
        self.thisptr = new pycolibri_classes.ClassEncoder(filename)

    def __dealloc__(self):
        del self.thisptr

    def __len__(self):
        return self.thisptr.size()

cdef class ClassDecoder:
    cdef pycolibri_classes.ClassDecoder *thisptr
    def __cinit__(self, str filename):
        self.thisptr = new pycolibri_classes.ClassDecoder(filename)

    def __dealloc__(self):
        del self.thisptr

    def __len__(self):
        return self.thisptr.size()

#cdef class Pattern:
#    cdef pycolibri.EncAnyGram * thisptr
#    def __init__(self):
#
#
#    def decode(ClassDecoder decoder):






