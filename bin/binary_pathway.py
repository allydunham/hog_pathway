#!/usr/bin/env python3
"""
Module for simulating biological signalling pathways using a binary on/off framework

ToDo
- Proteins to class?
- More sophisticated integration of upstream signal
- Add repression interactions?
- easier way to remove outputs
"""
class ProteinNetwork:
    """Signalling network"""
    def __init__(self, proteins=None, inputs=(), outputs=(),
                 threshold=0.5):
        self.proteins = {}
        if not proteins is None:
            for name, data in proteins.items():
                self.protein(name, **data)

        self.inputs = set(inputs)
        self.outputs = set(outputs)
        self.threshold = threshold

    def __repr__(self):
        return "ProteinNetwork(proteins={}, inputs={}, outputs={}, threshold={})\
                ".format(self.proteins, self.inputs, self.outputs, self.threshold)

    def __str__(self):
        # Potentially nicer representation of network
        return self.__repr__()

    def protein(self, name, outputs=None, inputs=None, function=None, overwrite=False):
        """Add a protein to the network or modify an exisiting one"""
        if not overwrite and name in self.proteins.keys():
            raise ValueError("Protein already exists and overwrite=False")
        elif overwrite and name in self.proteins.keys():
            if not outputs is None:
                self.proteins[name]['outputs'] = set(outputs)
            if not inputs is None:
                self.proteins[name]['inputs'] = set(inputs)
            if not function is None:
                self.proteins[name]['function'] = function
        else:
            if outputs is None:
                outputs = set()
            if inputs is None:
                inputs = set()
            if function is None:
                function = 1

            self.proteins[name] = {'inputs':set(inputs),
                                   'outputs':set(outputs),
                                   'function':function}
        self.verify_edges(name)

    def edge(self, source, sink):
        """Add an interaction between proteins already in the network"""
        self.proteins[source]['outputs'].append(sink)
        self.proteins[source]['inputs'].append(source)

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
        for i in self.proteins[name]['inputs']:
            try:
                self.proteins[i]['outputs'].add(name)
            except KeyError:
                print("Warning: absent protein {} referenced as input to {}".format(i, name))

        for i in self.proteins[name]['outputs']:
            try:
                self.proteins[i]['inputs'].add(name)
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
            if self.proteins[current_protein]['function'] > self.threshold:
                active[current_protein] = True
                to_test += list(self.proteins[current_protein]['outputs'])

        return any([active[x] for x in self.outputs])

if __name__ == "__main__":
    ## Can implement reduced version of hog with only direct interactions
    ## (i.e ignoring scaffold proteins and negative regulators)
    hog = ProteinNetwork()
    hog.protein('cdc42')
    hog.protein('sho1')
    hog.protein('sln1')

    hog.protein('cla4', inputs=['cdc42'])
    hog.protein('ste11', inputs=['cla4', 'sho1'])
    hog.protein('ypd1', inputs=['sln1'])
    hog.protein('ssk1', inputs=['ypd1'])
    hog.protein('ssk2', inputs=['ssk1'])
    hog.protein('pbs2', inputs=['ste11', 'ssk2'])
    hog.protein('hog1', inputs=['pbs2'])
    hog.protein('hot1', inputs=['hog1'])
    hog.protein('smp1', inputs=['hog1'])
    hog.protein('sko1', inputs=['hog1'])
    hog.protein('msn2', inputs=['hog1'])
    hog.protein('msn4', inputs=['hog1'])

    hog.set_input('cdc42', 'sho1', 'sln1')
    hog.set_output('hot1', 'smp1', 'msn2', 'msn4')
