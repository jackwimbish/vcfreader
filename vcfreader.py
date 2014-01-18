class VCFReader(object):
    """Class for opening and reading info from a VCF file."""

    def __init__(self, vcffile):
        """Constructor - read and process all lines up until the first data line.
        Accepts a file object to a vcf file or a file path.
        """
        if isinstance(vcffile, str):
            vcffile = open(vcffile,'r')
        self.filehandle = vcffile
        self.currentline = self.filehandle.readline()
        self.fileheader = ''
        while self.currentline[0:6] != '#CHROM':    # go to column header line
            self.fileheader += self.currentline
            self.currentline = self.filehandle.readline()
        # code below deals with column header line
        self.currentfields = self.currentline[1:-1].split('\t')
        if 'FORMAT' in self.currentfields:
            self.samples = self.currentfields[
                self.currentfields.index('FORMAT')+1:]
        self.columnindices = {}
        # map column names to indices
        for i in range(len(self.currentfields)):
            self.columnindices[self.currentfields[i]] = i
        self.setfields()
        self.last = None     # to prevent counting duplicate rows in the VCF
        self.unique = True   # this works only if vcf is coordinate sorted

    def getfieldval(self, fieldname):
        """Returns the value of a field with name fieldname."""
        return self.currentfields[self.columnindices[fieldname]]

    def setfields(self):
        """Set current fields to match current line of VCF file."""
        self.currentfields = self.currentline.strip().split('\t')
        self.CHROM = self.getfieldval('CHROM')
        self.POS = self.getfieldval('POS')
        self.ID = self.getfieldval('ID')
        self.REF = self.getfieldval('REF')
        self.ALT = self.getfieldval('ALT').split(',')
        self.QUAL = self.getfieldval('QUAL')
        self.FILTER = self.getfieldval('FILTER')
        if 'FORMAT' in self.columnindices:
            self.FORMAT = self.getfieldval('FORMAT').split(':')
        self.INFO = self.parseinfo()
        self.alleles = [self.REF]
        self.alleles.extend(self.ALT)

    def nextentry(self):
        """Move to next line of VCF file.  Returns that line."""
        self.currentline = self.filehandle.readline()
        if self.currentline:
            self.setfields()
            self.unique = self.checkunique()
        return self.currentline

    def checkunique(self):
        """Check that the current entry of the vcf file points to a
        new position.  Returns true for a new position, false for a repeat.
        This method requires that the vcf file is sorted by position to work
        correctly.
        """
        if self.last == (self.CHROM, self.POS, self.ALT):
            return False
        else:
            self.last = (self.CHROM, self.POS, self.ALT)
            return True

    def parseinfo(self):
        """Parse the INFO field of the current row of the vcf for data. Returns a
        dictionary containing keys with values for fields with a value. For
        boolean fields, the dictionary will have a value of True for the field
        name.
        """
        info_fields = self.getfieldval('INFO').split(';')
        info_map = {}
        for info_field in info_fields:
            parts = info_field.split('=')
            if len(parts) == 1:
                info_map[parts[0]] = True
            else:
                info_map[parts[0]] = parts[1]
        return info_map

    def getrgenotype(self, sampleid):
        """Get "relative genotype" of a sample where 0 is REF and 1 or greater
        are ALT geonotypes. (e.x. (0,1) for het)
        """
        sampledata = self.getsampinfo(sampleid)
        if not sampledata:
            return None
        rgenotype = sampledata[self.FORMAT.index('GT')]
        if '.' in rgenotype:
            return None
        delim = '/' if '/' in rgenotype else '|'
        return [int(x) for x in rgenotype.split(delim)]

    def getgenotype(self, sampleid):
        """Return the genotype of a sample."""
        rgenotype = self.getrgenotype(sampleid)
        if not rgenotype:
            return None
        return [self.alleles[x] for x in rgenotype]

    def getsampinfo(self, sampleid):
        """Return read information for a sample at the current position. Refer to
        the FORMAT string to interpret this information.
        """
        info = self.getfieldval(sampleid)
        if not info:
            return None
        if info == './.':
            None
        else:
            return info.split(':')

    def outputentry(self, newline=True):
        """Turns the current state of the vcfreader into a row for a vcf file.
        Returns row as a string.
        """
        cols = [self.CHROM, self.POS, self.ID, self.REF, ','.join(self.ALT),
                self.QUAL, self.FILTER]
        info_fields = []
        for key, val in sorted(self.INFO.items()):
            if val is True:
                info_fields.append(key)
            else:
                info_fields.append('{0}={1}'.format(key, val))
        cols.append(';'.join(info_fields))
        cols.append(':'.join(self.FORMAT))
        cols.extend(self.getfieldval(samp) for samp in self.samples)
        cols = [str(col) for col in cols]
        return '\t'.join(cols) + ('\n' if newline else '')

    def colheaders(self, newline=True):
        """Generates a column header line as a string."""
        cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
                'FORMAT']
        cols.extend(self.samples)
        return '#' + '\t'.join(cols) + ('\n' if newline else '')

    def reset(self):
        """Restarts at beginning of VCF file."""
        self.filehandle.seek(0)
        self.__init__(self.filehandle)

