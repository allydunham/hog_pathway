#!/usr/bin/env python3
"""
Module for simulating biological signalling pathways using a binary on/off framework

ToDo
- More sophisticated integration of upstream signal
- Add repression interactions?
- easier way to remove outputs
- deal with loop interactions
"""
FUNCTIONALITY_THRESHOLD = 0.5

class ProteinNetwork:
    """Signalling network"""
    def __init__(self, inputs=(), outputs=()): #, proteins=None
        self.proteins = {}
        # if not proteins is None:
        #     for name, data in proteins.items():
        #         self.node(name, **data)

        self.inputs = set(inputs)
        self.outputs = set(outputs)

    def __repr__(self):
        return "ProteinNetwork(proteins={}, inputs={}, outputs={}, threshold={})\
                ".format(self.proteins, self.inputs, self.outputs)

    def __str__(self):
        # Potentially nicer representation of network
        return self.__repr__()

    def node(self, obj):
        """Add a protein to the network or modify an exisiting one"""
        self.proteins[obj.name] = obj
        self.verify_edges(obj.name)

    def edge(self, source, sink):
        """Add an interaction between proteins already in the network"""
        self.proteins[source].add_output(sink)
        self.proteins[sink].add_input(source)

    def set_input(self, *names):
        """Set a protein(s) as an input source"""
        names = set(names)
        if all([x in self.proteins.keys() for x in names]):
            self.inputs.update(names)
        else:
            raise ValueError("Protein(s) not in network")

    def set_output(self, *names):
        """Set a protein(s) as an output source"""
        names = set(names)
        if all([x in self.proteins.keys() for x in names]):
            self.outputs.update(names)
        else:
            raise ValueError("Protein(s) not in network")

    def verify_edges(self, name):
        """Check edges are fully described with respect to named protein"""
        for i in self.proteins[name].inputs:
            try:
                self.proteins[i].add_output(name)
            except KeyError:
                print("Warning: absent protein {} referenced as input to {}".format(i, name))

        for i in self.proteins[name].outputs:
            try:
                self.proteins[i].add_input(name)
            except KeyError:
                print("Warning: absent protein {} referenced as output to {}".format(i, name))

    def verify_all_edges(self):
        """Check edges are fully described across all proteins"""
        for protein in self.proteins:
            self.verify_edges(protein)

    def get_activity(self):
        """Check activity of the network"""
        active = {name: False for name in self.proteins}
        to_test = list(self.inputs)
        while to_test:
            current_protein = to_test.pop()
            if self.proteins[current_protein].get_activity():
                active[current_protein] = True
                to_test += list(self.proteins[current_protein].outputs)

        return any([active[x] for x in self.outputs])

class Protein:
    """Individual component of a protein network"""
    def __init__(self, name, inputs=None, outputs=None, function=1):
        if inputs is None:
            inputs = set()
        if outputs is None:
            outputs = set()

        self.name = name
        self.inputs = set(inputs)
        self.outputs = set(outputs)
        self.function = function

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    def add_input(self, name):
        """Add an input edge"""
        self.inputs.add(name)

    def add_output(self, name):
        """Add an output edge"""
        self.outputs.add(name)

    def get_activity(self):
        """Test if the protein is functional"""
        return self.function > FUNCTIONALITY_THRESHOLD

class ProteinComplex(Protein):
    """Protein network component made up of multiple individual components"""
    def __init__(self, name, components, inputs=None, outputs=None):
        Protein.__init__(self, name, inputs, outputs, None)
        self.components = {name: Protein(name, inputs=None, outputs=None, function=func)
                           for name, func in components.items()}

    def get_activity(self):
        """Test if the complex is active"""
        return all([x.get_activity() for x in self.components.values()])


if __name__ == "__main__":
    ## Can implement reduced version of hog with only direct interactions
    ## (i.e ignoring scaffold proteins and negative regulators)
    hog = ProteinNetwork()
    hog.node(Protein('cdc42'))
    hog.node(ProteinComplex('sho1-hkr1-msb1', {'sho1': 1, 'hkr1': 1, 'msb1': 1}))
    hog.node(Protein('sln1'))

    hog.node(ProteinComplex('cla4-ste20', {'cla4':1, 'ste20': 1}, inputs=['cdc42']))
    hog.node(Protein('ste11', inputs=['cla4-ste20', 'sho1-hkr1-msb1']))
    hog.node(Protein('ypd1', inputs=['sln1']))
    hog.node(Protein('ssk1', inputs=['ypd1']))
    hog.node(ProteinComplex('ssk2-ssk22', {'ssk2': 1, 'ssk22': 1}, inputs=['ssk1']))
    hog.node(Protein('pbs2', inputs=['ste11', 'ssk2-ssk22']))
    hog.node(Protein('hog1', inputs=['pbs2']))
    hog.node(Protein('hot1', inputs=['hog1']))
    hog.node(Protein('smp1', inputs=['hog1']))
    hog.node(Protein('sko1', inputs=['hog1']))
    hog.node(Protein('msn2', inputs=['hog1']))
    hog.node(Protein('msn4', inputs=['hog1']))

    hog.set_input('cdc42', 'sho1-hkr1-msb1', 'sln1')
    hog.set_output('hot1', 'smp1', 'msn2', 'msn4')
