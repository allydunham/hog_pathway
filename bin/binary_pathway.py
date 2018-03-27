#!/usr/bin/env python3
"""
Module for simulating biological signalling pathways using a binary on/off framework

ToDo
- More sophisticated integration of upstream signal
- Add repression interactions?
- easier way to remove outputs
- deal with loop interactions
- More complex functionality representation
"""
FUNCTIONALITY_THRESHOLD = 0.05

class ProteinNetwork:
    """Signalling network"""
    def __init__(self, proteins=None, inputs=(), outputs=()):
        self.nodes = {}
        self.complex_members = {} # Lookup of all protein:complex pairs
        if not proteins is None:
            for protein in proteins.values():
                self.node(protein)
                if isinstance(protein, ProteinComplex):
                    for name in protein.components:
                        self.complex_members[name] = protein.name

        self.inputs = set(inputs)
        self.outputs = set(outputs)


    def __repr__(self):
        return "ProteinNetwork(proteins={}, inputs={}, outputs={})".format(
            self.nodes, self.inputs, self.outputs
            )

    def __str__(self):
        # Potentially nicer representation of network
        return self.__repr__()

    def node(self, obj):
        """Add a protein to the network or modify an exisiting one"""
        self.nodes[obj.name] = obj

        if isinstance(obj, ProteinComplex):
            for name in obj.components:
                self.complex_members[name] = obj.name

        self.verify_edges(obj.name)

    def edge(self, source, sink):
        """Add an interaction between proteins already in the network"""
        self.nodes[source].add_output(sink)
        self.nodes[sink].add_input(source)

    def set_input(self, *names, reset=False):
        """Set a protein(s) as an input source"""
        names = set(names)
        if not names or all([x in self.nodes.keys() for x in names]):
            if reset:
                self.inputs = names
            else:
                self.inputs.update(names)
        else:
            raise ValueError("Protein(s) not in network")

    def set_output(self, *names, reset=False):
        """Set a protein(s) as an output source"""
        names = set(names)
        if not names or all([x in self.nodes.keys() for x in names]):
            if reset:
                self.outputs = names
            else:
                self.outputs.update(names)
        else:
            raise ValueError("Protein(s) not in network")

    def set_function(self, name, function):
        """Set a component proteins functionality"""
        if name in self.nodes:
            self.nodes[name].function = function
        elif name in self.complex_members:
            comp = self.complex_members[name]
            self.nodes[comp].components[name].function = function
        else:
            raise ValueError("Protein not in network")

    def verify_edges(self, name):
        """Check edges are fully described with respect to named protein"""
        for i in self.nodes[name].inputs:
            try:
                self.nodes[i].add_output(name)
            except KeyError:
                print("Warning: absent protein {} referenced as input to {}".format(i, name))

        for i in self.nodes[name].outputs:
            try:
                self.nodes[i].add_input(name)
            except KeyError:
                print("Warning: absent protein {} referenced as output to {}".format(i, name))

    def verify_all_edges(self):
        """Check edges are fully described across all proteins"""
        for name in self.nodes:
            self.verify_edges(name)

    def get_activity(self):
        """Check activity of the network"""
        active = {name: False for name in self.nodes}
        to_test = list(self.inputs)
        while to_test:
            current_protein = to_test.pop()
            if self.nodes[current_protein].get_activity():
                active[current_protein] = True
                to_test += list(self.nodes[current_protein].outputs)

        return any([active[x] for x in self.outputs])

    def get_protein_names(self):
        """Return a list of all proteins in the network"""
        proteins = []
        for name, obj in self.nodes.items():
            if isinstance(obj, ProteinComplex):
                proteins.extend(obj.components.keys())
            elif isinstance(obj, Protein):
                proteins.append(name)
            else:
                # Should never occur
                raise ValueError("Non-protein/complex object in network")

        return proteins

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
        return "Protein('{}', inputs={}, outputs={}, function={})".format(
            self.name, self.inputs, self.outputs, self.function
            )

    def __str__(self):
        return self.__repr__()

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
        self.components = {x.name: x for x in components}

    def __repr__(self):
        return "ProteinComplex('{}', {},inputs={}, outputs={})".format(
            self.name, [x for x in self.components.values()], self.inputs, self.outputs
            )

    def get_activity(self):
        """Test if the complex is active"""
        return all([x.get_activity() for x in self.components.values()])


if __name__ == "__main__":
    ## Can implement reduced version of hog without negative regulators (e.g. phosphatases)
    hog = ProteinNetwork()
    hog.node(Protein('cdc42'))
    hog.node(ProteinComplex('sho1-hkr1-msb1',
                            [Protein('sho1'), Protein('hkr1'), Protein('msb1')]))

    hog.node(Protein('sln1'))
    hog.node(ProteinComplex('cla4-ste20',
                            [Protein('cla4'), Protein('ste20')],
                            inputs=['cdc42']))

    hog.node(Protein('ste11', inputs=['cla4-ste20', 'sho1-hkr1-msb1']))
    hog.node(Protein('ypd1', inputs=['sln1']))
    hog.node(Protein('ssk1', inputs=['ypd1']))
    hog.node(ProteinComplex('ssk2-ssk22',
                            [Protein('ssk2'), Protein('ssk22')],
                            inputs=['ssk1']))

    hog.node(Protein('pbs2', inputs=['ste11', 'ssk2-ssk22']))
    hog.node(Protein('hog1', inputs=['pbs2']))
    hog.node(Protein('hot1', inputs=['hog1']))
    hog.node(Protein('smp1', inputs=['hog1']))
    hog.node(Protein('sko1', inputs=['hog1']))
    hog.node(Protein('msn2', inputs=['hog1']))
    hog.node(Protein('msn4', inputs=['hog1']))

    hog.set_input('cdc42', 'sho1-hkr1-msb1', 'sln1')
    hog.set_output('hot1', 'smp1', 'msn2', 'msn4')
