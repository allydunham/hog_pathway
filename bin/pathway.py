#!/usr/bin/env python3
"""
Module for simulating biological signalling pathways using a binary on/off framework

ToDo
- More sophisticated integration of upstream signal
- Add repression interactions?
- easier way to remove outputs
- deal with loop interactions
- More complex functionality representation
- Separate binding/localisation and kinase interactions
"""
import numpy as np
FUNCTIONALITY_THRESHOLD = 0.5

class ProteinNetwork:
    """Signalling network"""
    def __init__(self, proteins=None, inputs=(), outputs=()):
        self.nodes = {} # dict of named nodes in network
        self.complex_members = {} # Lookup of all protein:complex pairs

        # Import any supplied nodes
        if not proteins is None:
            for protein in proteins.values():
                self.node(protein)
                if isinstance(protein, ProteinComplex):
                    for name in protein.components:
                        self.complex_members[name] = protein.name

        self.inputs = set(inputs) # Set of input nodes (activity sources)
        self.outputs = set(outputs) # Set of output nodes (netowrk targets)
        self.branches = [] # list of branches in the network

        self.verify_all_edges()

    def __repr__(self):
        return "ProteinNetwork(proteins={}, inputs={}, outputs={})".format(
            self.nodes, self.inputs, self.outputs
            )

    def __str__(self):
        # Potentially nicer representation of network
        return self.__repr__()

    def node(self, obj):
        """Add a node to the network or modify an exisiting one"""
        self.nodes[obj.name] = obj

        if isinstance(obj, ProteinComplex):
            for name in obj.components:
                self.complex_members[name] = obj.name

        self.verify_edges(obj.name)

    def edge(self, source, sink):
        """Add an interaction between nodes already in the network"""
        self.nodes[source].add_output(sink)
        self.nodes[sink].add_input(source)

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

    def set_probability(self, name, probability):
        """Set a component proteins functionality probability"""
        if name in self.complex_members: # Set complex first to avoid setting directly
            comp = self.complex_members[name]
            self.nodes[comp].components[name].probability_functional = probability
        elif name in self.nodes:
            self.nodes[name].probability_functional = probability
        else:
            raise ValueError("Protein not in network")

    def calc_activity_probabilities(self):
        """Determine the probability that each network component is active based\
           on upstream activity probabilities"""
        # Reset probabilities
        for node in self.nodes.values():
            node.probability_active = np.nan

        # Set input probabilities and initial downstream nodes
        stack = []
        for name in self.inputs:
            node = self.nodes[name]
            node.probability_active = node.get_probability_functional()
            stack.extend(node.outputs)

        # Iterate over nodes to get probability
        while stack:
            name = stack.pop()
            node = self.nodes[name]
            upstream_p = [self.nodes[n].probability_active for n in
                          node.inputs]

            # if any upstream prob isn't available move node to bottom of stack
            if any(np.isnan(p) for p in upstream_p):
                stack.insert(0, name)
                continue

            else:
                # Calc active p as p_func * (1 - prod(not outputs active))
                node.probability_active = (node.get_probability_functional() *
                                           (1 - np.prod([1 - x for x in upstream_p])))

                # Add outputs to stack
                for i in node.outputs:
                    if not i in stack:
                        stack.insert(0, i)

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

    def get_activity(self, threshold=FUNCTIONALITY_THRESHOLD):
        """Check binary activity of the network outputs based on P(Active) > threshold"""
        active = {name: False for name in self.nodes}
        to_test = list(self.inputs)
        while to_test:
            current_protein = to_test.pop()
            if self.nodes[current_protein].get_activity(threshold):
                active[current_protein] = True
                to_test += list(self.nodes[current_protein].outputs)

        return any([active[x] for x in self.outputs])

class Protein:
    """Individual component of a protein network"""
    def __init__(self, name, inputs=None, outputs=None, p_functional=0.95):
        if inputs is None:
            inputs = set()
        if outputs is None:
            outputs = set()

        self.name = name
        self.inputs = set(inputs)
        self.outputs = set(outputs)
        self.probability_functional = p_functional
        self.probability_active = np.nan


    def __repr__(self):
        return "Protein('{}', inputs={}, outputs={}, p_functional={})".format(
            self.name, self.inputs, self.outputs, self.probability_functional)

    def __str__(self):
        return self.__repr__()

    def add_input(self, name):
        """Add an input edge"""
        self.inputs.add(name)

    def add_output(self, name):
        """Add an output edge"""
        self.outputs.add(name)

    def get_activity(self, threshold=FUNCTIONALITY_THRESHOLD):
        """Test if the protein is functional"""
        return self.probability_functional > threshold

    def get_probability_functional(self):
        """Get the probability of the protein being functional"""
        return self.probability_functional

class ProteinComplex(Protein):
    """Protein network component made up of multiple individual components"""
    def __init__(self, name, components, inputs=None, outputs=None):
        Protein.__init__(self, name, inputs, outputs, None)
        self.components = {x.name: x for x in components}

    def __repr__(self):
        return "ProteinComplex('{}', {}, inputs={}, outputs={})".format(
            self.name, [x for x in self.components.values()], self.inputs, self.outputs)

    def get_activity(self, threshold=FUNCTIONALITY_THRESHOLD):
        """Test if the complex is active"""
        return all([x.get_activity(threshold) for x in self.components.values()])

    def get_probability_functional(self):
        """Get the probability of the complex being functional"""
        return np.prod([x.get_probability_functional() for x in self.components.values()])


if __name__ == "__main__":
    ## Can implement reduced version of hog without negative regulators (e.g. phosphatases)
    hog = ProteinNetwork()
    hog.node(ProteinComplex('sho1-hkr1-msb1',
                            [Protein('sho1'), Protein('hkr1'), Protein('msb1')]))

    hog.node(Protein('sln1'))

    hog.node(ProteinComplex('cla4-ste20-cdc42',
                            [Protein('cla4'), Protein('ste20'), Protein('cdc42')],
                            inputs=['sho1-hkr1-msb1']))

    hog.node(ProteinComplex('ste11-ste50',
                            [Protein('ste11'), Protein('ste50')],
                            inputs=['cla4-ste20-cdc42']))

    hog.node(Protein('ypd1', inputs=['sln1']))
    hog.node(Protein('ssk1', inputs=['ypd1']))
    hog.node(ProteinComplex('ssk2-ssk22',
                            [Protein('ssk2'), Protein('ssk22')],
                            inputs=['ssk1']))

    hog.node(Protein('pbs2', inputs=['ste11-ste50', 'ssk2-ssk22']))
    hog.node(Protein('hog1', inputs=['pbs2']))
    hog.node(Protein('hot1', inputs=['hog1']))
    hog.node(Protein('smp1', inputs=['hog1']))
    hog.node(Protein('sko1', inputs=['hog1']))
    hog.node(Protein('msn2', inputs=['hog1']))
    hog.node(Protein('msn4', inputs=['hog1']))

    hog.set_input('sho1-hkr1-msb1', 'sln1')
    hog.set_output('hot1', 'smp1', 'msn2', 'msn4')

    print(hog)
    hog.calc_activity_probabilities()
