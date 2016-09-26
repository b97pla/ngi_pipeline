class NGIObject(object):
    def __init__(self, name, dirname, subitem_type):
        self.name = name
        self.dirname = dirname
        self._subitems = {}
        self._subitem_type = subitem_type

    def _add_subitem(self, name, dirname):
        # Only add a new item if the same item doesn't already exist
        try:
            subitem = self._subitems[name]
        except KeyError:
            subitem = self._subitems[name] = self._subitem_type(name, dirname)
        return subitem

    def __iter__(self):
        return iter(self._subitems.values())

    def __unicode__(self):
        return self.name

    def __str__(self):
        return self.__unicode__()

    def __repr__(self):
        return "{}: \"{}\"".format(type(self), self.name)


## TODO consider changing the default __repr__ and __str__ to project_id
class NGIProject(NGIObject):
    def __init__(self, name, dirname, project_id, base_path):
        self.base_path = base_path
        super(NGIProject, self).__init__(name, dirname, subitem_type=NGISample)
        self.samples = self._subitems
        self.add_sample = self._add_subitem
        self.project_id = project_id
        self.command_lines = []
        self.chip_genotypes = None


class NGISample(NGIObject):
    def __init__(self, *args, **kwargs):
        super(NGISample, self).__init__(subitem_type=NGILibraryPrep, *args, **kwargs)
        self.libpreps = self._subitems
        self.add_libprep = self._add_subitem


class NGILibraryPrep(NGIObject):
    def __init__(self, *args, **kwargs):
        super(NGILibraryPrep, self).__init__(subitem_type=NGISeqRun, *args, **kwargs)
        self.seqruns = self._subitems
        self.add_seqrun = self._add_subitem


class NGISeqRun(NGIObject):
    def __init__(self, *args, **kwargs):
        super(NGISeqRun, self).__init__(subitem_type=None, *args, **kwargs)
        self.fastq_files = self._subitems = []
        ## Not working
        #delattr(self, "_add_subitem")

    def __iter__(self):
        return iter(self._subitems)

    def add_fastq_files(self, fastq):
        if type(fastq) == list:
            self._subitems.extend(fastq)
        elif type(fastq) == str or type(fastq) == unicode:
            self._subitems.append(str(fastq))
        else:
            raise TypeError("Fastq files must be passed as a list or a string: " \
                            "got \"{}\"".format(fastq))

class NGIChipGenotypes(NGIObject):
    def __init__(self, *args, **kwargs):
        super(NGIChipGenotypes, self).__init__(*args, **kwargs, subitem_type=None)

