import msprime
import tskit

# NOTE: the Python edge differences API
# is pretty feature rich.
# In lower level languages, some operations
# require some extra manual book-keeping


def mutations_in_each_tree(ts: tskit.TreeSequence):
    current_site_index = 0
    current_tree = 0
    for e in ts.edge_diffs():
        while (
            current_site_index < ts.num_sites
            and ts.site(current_site_index).position < e.interval.right
        ):
            print(f"site {ts.site(current_site_index)} is in tree {current_tree}")
            current_site_index += 1
        current_tree += 1
    assert current_tree == ts.num_trees


def parents_in_each_tree(ts: tskit.TreeSequence):
    parents = [tskit.NULL for _ in range(ts.num_nodes + 1)]
    for e, t in zip(ts.edge_diffs(), ts.trees()):
        for o in e.edges_out:
            parents[o.child] = tskit.NULL
        for i in e.edges_in:
            parents[i.child] = i.parent
        assert all([i == j for i, j in zip(parents, t.parent_array)])


def num_children_of_each_parent(ts: tskit.TreeSequence):
    num_children = [0] * ts.num_nodes
    for e, t in zip(ts.edge_diffs(), ts.trees()):
        for o in e.edges_out:
            num_children[o.parent] -= 1
            assert num_children[o.parent] >= 0
        for i in e.edges_in:
            num_children[i.parent] += 1
        for i, j in enumerate(num_children):
            assert t.num_children(i) == j


def num_samples_descending_from_each_node(ts: tskit.TreeSequence):
    # NOTE: once we know the number of samples descending from a node
    # in a given tree AND we know how to iterate over all mutations
    # in a tree, we can ask how many samples are below a mutation!!!!
    parents = [tskit.NULL for _ in range(ts.num_nodes + 1)]
    num_samples = [0] * ts.num_nodes
    for s in ts.samples():
        num_samples[s] += 1

    for e, t in zip(ts.edge_diffs(), ts.trees(sample_lists=True)):
        for o in e.edges_out:
            p = parents[o.child]
            while p != tskit.NULL:
                num_samples[p] -= num_samples[o.child]
                p = parents[p]
            parents[o.child] = tskit.NULL
        for i in e.edges_in:
            p = i.parent
            while p != tskit.NULL:
                num_samples[p] += num_samples[i.child]
                p = parents[p]
            parents[i.child] = i.parent
        for i, j in enumerate(num_samples):
            assert t.num_samples(i) == j, f"{t.num_samples(i)} != {j}"


ts = msprime.sim_ancestry(
    5,
    population_size=10000,
    sequence_length=1e8,
    recombination_rate=1e-10,
    random_seed=666,
)
assert ts.num_trees > 0
print(ts.num_trees)

ts = msprime.sim_mutations(ts, rate=1e-10, random_seed=666)
assert ts.num_mutations > 0
print(ts.num_mutations)

mutations_in_each_tree(ts)
parents_in_each_tree(ts)
num_children_of_each_parent(ts)
num_samples_descending_from_each_node(ts)
