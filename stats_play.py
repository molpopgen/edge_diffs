import tskit

tables = tskit.TableCollection(10.0)

n0 = tables.nodes.add_row(0, time=2)
n1 = tables.nodes.add_row(0, time=1)
n2 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0)
n3 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0)
n4 = tables.nodes.add_row(tskit.NODE_IS_SAMPLE, time=0)

tables.edges.add_row(0.0, tables.sequence_length, n1, n2)
tables.edges.add_row(0.0, tables.sequence_length, n1, n3)
tables.edges.add_row(0.0, tables.sequence_length, n0, n1)
tables.edges.add_row(0.0, tables.sequence_length, n0, n4)


site = tables.sites.add_row(5, ancestral_state="A")
tables.mutations.add_row(site, node=n1, time=1.0, derived_state="G")
tables.mutations.add_row(site, node=n3, time=0.1, derived_state="C")
tables.mutations.add_row(site, node=n1, time=1.1, derived_state="T")

tables.sort()

ts = tables.tree_sequence()

print(ts.draw_text())

print(ts.tables.mutations)

print(ts.diversity(span_normalise=False))

for g in ts.haplotypes():
    print(g)

tables.compute_mutation_parents()
ts = tables.tree_sequence()
print(ts.diversity(span_normalise=False))

parent = [tskit.NULL] * ts.num_nodes
num_samples_below = [tskit.NULL] * ts.num_nodes
for diffs in ts.edge_diffs():
    print(diffs)
